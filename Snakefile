configfile: "config.yaml"

# Liste des échantillons (clés de config["samples"] : SRR -> condition)
SAMPLES = list(config["samples"])

# Index Bowtie : bowtie-build produit 6 fichiers, Snakemake doit tous les connaître
INDEX_PREFIX = "results/index/reference"
INDEX_FILES = multiext(
    INDEX_PREFIX,
    ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt",
    ".rev.1.ebwt", ".rev.2.ebwt",
)

# Cible finale : les deux figures à reproduire
#   - Supp. Fig. 3 : MA-plot de tout le jeu de données
#   - Fig. 3c      : MA-plot des gènes liés à la traduction
rule all:
    input:
        "results/figures/MA_plot_all.png",
        "results/figures/MA_plot_translation.png"

# Cible du rendu intermédiaire (13 novembre) : données + génome + index
#   lancer avec :  snakemake --use-apptainer --cores 4 midterm
rule midterm:
    input:
        expand("results/trimmed/{sample}.trimmed.fastq.gz", sample=SAMPLES),
        INDEX_FILES,
        "data/reference.gff"
        

rule download:
    """Télécharge les FASTQ depuis la SRA (reads single-end) puis les compresse."""
    output: "data/{sample}.fastq.gz"
    threads: 4
    container: "containers/sratoolkit.sif"
    shell:
        """
        fasterq-dump --threads {threads} --progress {wildcards.sample} -O data/
        gzip data/{wildcards.sample}.fastq
        """


rule trim:
    input: "data/{sample}.fastq.gz"
    output: "results/trimmed/{sample}.trimmed.fastq.gz"
    container: "containers/cutadapt.sif"
    shell: "cutadapt -q 20 -m 25 -o {output} {input}"


rule get_genome:
    """Télécharge le génome de référence (S. aureus NCTC 8325, CP000253.1)."""
    output: "data/reference.fasta"
    container: "containers/utils.sif" 
    shell:
        """
        wget -q -O {output} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        """

rule get_annotation:
    """Télécharge l'annotation GFF3 du même génome."""
    output: "data/reference.gff"
    container: "containers/utils.sif"
    shell:
        """
        wget -q -O {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
        """


rule bowtie_index:
    """Construit l'index Bowtie à partir du fasta (une seule fois, partagé par tous les échantillons)."""
    input: "data/reference.fasta"
    output: INDEX_FILES
    params: prefix=INDEX_PREFIX   # préfixe (texte), pas un fichier : bowtie-build crée les 6 fichiers
    log: "logs/bowtie_index.log"
    container: "containers/bowtie_1.3.1.img"   # écart avec l'article (0.12.7), à justifier dans le rapport
    shell: "bowtie-build-s {input} {params.prefix} > {log} 2>&1"


rule bowtie_map:
    """Mappe les lectures trimmées avec Bowtie (sortie SAM)."""
    input:
        reads="results/trimmed/{sample}.trimmed.fastq.gz",
        index=INDEX_FILES   # le mapping ne démarre qu'une fois l'index construit
    output: temp("results/mapping/{sample}.sam")   # temp : supprimé quand plus nécessaire
    params: prefix=INDEX_PREFIX
    log: "logs/bowtie_map_{sample}.log"   # Bowtie y écrit son taux de mapping
    threads: 4
    container: "containers/bowtie_1.3.1.img"
    shell:
        "bowtie-align-s -p {threads} -S {params.prefix} <(gunzip -c {input.reads}) > {output} 2> {log}"


rule samtools_sort:
    """Trie le SAM, produit le BAM, puis l'indexe."""
    input: "results/mapping/{sample}.sam"
    output:
        bam="results/mapping/{sample}.bam",
        bai="results/mapping/{sample}.bam.bai"
    threads: 4
    container: "containers/samtools_1.20.0.img"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        """

rule count:
    """featureCounts 1.4.6-p3 : comptage par gène, library stranded (-s 1)."""
    input:
        bams=expand("results/mapping/{sample}.bam", sample=SAMPLES),
        gff="data/reference.gff"
    output: "results/counts/counts.txt"
    threads: 4
    log: "logs/featurecounts.log"
    container: "containers/subread.sif"
    shell:
        """
        featureCounts --extraAttributes Name -t gene -g ID -s 1 -F GTF \
            -T {threads} -a {input.gff} -o {output} {input.bams} > {log} 2>&1
        """

rule translation_genes:
    """Liste des gènes de la traduction via l'API REST de KEGG (voir config.yaml)."""
    output: "results/annotation/translation_genes.txt"
    params: pathways=config["kegg_pathways"]
    container: "containers/utils.sif"
    shell: "python3 scripts/kegg_translation_genes.py {params.pathways} > {output}"

rule deseq2:
    """DESeq2 (persister vs control) puis MA-plots (tout le jeu / gènes de la traduction)."""
    input:
        counts="results/counts/counts.txt",
        translation="results/annotation/translation_genes.txt"
    output:
        table="results/deseq2/deseq2_results.csv",
        ma_all="results/figures/MA_plot_all.png",
        ma_translation="results/figures/MA_plot_translation.png"
    params:
        # "SRR1:persister,SRR2:control,..." construit depuis config.yaml
        samples=",".join(f"{s}:{c}" for s, c in config["samples"].items())
    log: "logs/deseq2.log"
    container: "containers/deseq2.sif"
    shell:
        """
        Rscript scripts/deseq2_ma.R {input.counts} {input.translation} \
            {params.samples} results > {log} 2>&1
        """