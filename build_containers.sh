#!/bin/bash
# Prépare tous les conteneurs du projet.
set -e
mkdir -p containers

# 1. Images téléchargées depuis Damona (Zenodo)
wget -nc -O containers/samtools_1.20.0.img \
  https://zenodo.org/record/23326397/files/samtools_1.20.0.img

# 2. Images construites depuis nos recettes .def (même nom pour .def et .sif)
for name in cutadapt sratoolkit utils subread deseq2; do
    if [ ! -f "containers/${name}.sif" ]; then
        apptainer build "containers/${name}.sif" "containers/${name}.def"
    fi
done

# 3. Bowtie : le .sif porte la version dans son nom, donc hors de la boucle
[ -f containers/bowtie_0.12.7.sif ] || apptainer build containers/bowtie_0.12.7.sif containers/bowtie.def