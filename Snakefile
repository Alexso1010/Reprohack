configfile: "config.yaml"

# Index Bowtie : on doit générer 6 fichiers, snakemake doit tous les connaitre
INDEX_PREFIX = "results/index/reference"
INDEX_FILES = multiext(
	INDEX_PREFIX,
	".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt",
	".rev.1.ebwt", ".rev.2.ebwt",
)

Rule all:
	input:
		expand("results/mapping/{sample}.bam.bai, sample=config["samples"])
# J'ai retiré la règle all de Thomas par ce que celle du dessus appelle la même chose que l'ancienne à la suite

rule download:
    output: "data/{sample}.fastq.gz"
    threads: 4
    container: "containers/sratoolskit.sif"
    shell: """
	fasterq-dump --threads {threads} --progress {wildcards.sample} -O data/
        gzip data/{wildcards.sample}.fastq
        """

rule trim:
    input: "data/{sample}.fastq.gz"
    output: "results/trimmed/{sample}.trimmed.fastq.gz"
    container: "containers/cutadapt.sif"
    shell: "cutadapt -q 20 -m 25 -o {output} {input}


# Télécharger le génome de référence
rule get_genome:
	output: "data/reference.fasta"
	container: "containers/bowtie.sif" #vérifier si wget dans le container concerné
	shell:
		"""
        	wget -q -O {output} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        	"""


rule bowtie_index:
	"""construction de l'index bowtie à partie du fasta"""
	input= "data/reference.fasta"
	output: INDEX_FILES
	params: prefix=INDEX_PREFIX #valeur fixe passée à la commande (ni fichier d'entrée ni de sortie) pour que la règle retrouve les fichiers elle même 
	log: "logs/bowtie_index.log"
	container: "container/bowtie_1.3.1.img" #version de bowtie à modifier
	shell: "bowtie-build {input} {params.prefix} > {log} 2>&1"
# Pas de {sample}, l'index est construit une seule fois seulement puis partagé par tous les échantillons


rule bowtie_map:
	"""Mapping des lectures triées, tri, et indexation du BAM"
	input: 
		reads="results/trimmed/{sample}.trimmed.fastq.gz", #appel des entrées nommés dans le shell
		index=INDEX_FILES #démarrage du mapping seulement une fois que l'index est construit
	output:
		bam="results/mapping/{sample}.bam"
		bai="results/mapping/{sample}.bam.bai"
	params : prefix=INDEX_PREFIX
	log: "logs/bowtie_map_{sample}.log" 	#conservation du taux de mapping 
	threads: 4	#nb de coeurs réservés
	container: "containers/bowtie.sif" #vérifier que container contient bowtie et samtools
	shell:
		"""
		bowtie -p {threads} -S {params.prefix} <(gunzip -c {input.reads}) 2> {log} \
            		| samtools sort -@ {threads} > {output.bam}
        	samtools index {output.bam}
		"""
	
