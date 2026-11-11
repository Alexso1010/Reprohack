configfile: "config.yaml"

# Index Bowtie : bowtie-build produit 6 fichiers, Snakemake doit tous les connaître
INDEX_PREFIX = "results/index/reference"
INDEX_FILES = multiext(
    INDEX_PREFIX,
    ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt",
    ".rev.1.ebwt", ".rev.2.ebwt",
)

# Cible finale : demander les .bai déclenche toute la chaîne
# (download -> trim -> mapping -> tri/index)
rule all:
    input:
        expand("results/mapping/{sample}.bam.bai", sample=config["samples"])


rule download:
    output: "data/{sample}.fastq.gz"
    threads: 4
    container: "containers/sratoolskit.sif"
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
    container: "containers/bowtie_1.3.1.img"   # À VÉRIFIER : wget présent dans l'image ?
    shell:
        """
        wget -q -O {output} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        """


rule bowtie_index:
    """Construit l'index Bowtie à partir du fasta (une seule fois, partagé par tous les échantillons)."""
    input: "data/reference.fasta"
    output: INDEX_FILES
    params: prefix=INDEX_PREFIX   # préfixe (texte), pas un fichier : bowtie-build crée les 6 fichiers
    log: "logs/bowtie_index.log"
    container: "containers/bowtie_1.3.1.img"   # écart avec l'article (0.12.7), à justifier dans le rapport
    shell: "bowtie-build {input} {params.prefix} > {log} 2>&1"


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
        "bowtie -p {threads} -S {params.prefix} <(gunzip -c {input.reads}) > {output} 2> {log}"


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
