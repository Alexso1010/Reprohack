#!/bin/bash
# Prépare tous les conteneurs du projet.
set -e
mkdir -p containers

# 1. Images téléchargées depuis Damona (Zenodo)
wget -nc -O containers/bowtie_1.3.1.img \
  https://zenodo.org/record/5810335/files/bowtie_1.3.1.img
wget -nc -O containers/samtools_1.20.0.img \
  https://zenodo.org/record/23326397/files/samtools_1.20.0.img

# 2. Images construites depuis nos recettes .def
for name in cutadapt sratoolkit utils subread deseq2; do
    if [ ! -f "containers/${name}.sif" ]; then
        apptainer build "containers/${name}.sif" "containers/${name}.def"
    fi
done