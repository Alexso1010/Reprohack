configfile: "config.yaml"

rule all:
    input:
        expand("results/trimmed/{sample}.trimmed.fastq.gz", sample=config["samples"])

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
    shell: "cutadapt -q 20 -m 25 -o {output} {input}"
