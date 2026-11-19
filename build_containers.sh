#!/bin/bash
# Prépare tous les conteneurs du projet.
set -e
mkdir -p containers

# 2. Images construites depuis nos recettes .def (même nom pour .def et .sif)
for name in cutadapt sratoolkit utils subread deseq2 mapping; do
    if [ ! -f "containers/${name}.sif" ]; then
        apptainer build "containers/${name}.sif" "containers/${name}.def"
    fi
done