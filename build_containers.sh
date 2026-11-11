#!/bin/bash
# Télécharge les images Apptainer depuis Damona (Zenodo)
set -e
mkdir -p containers

wget -nc -O containers/bowtie_1.3.1.img \
  https://zenodo.org/record/5810335/files/bowtie_1.3.1.img

wget -nc -O containers/samtools_1.20.0.img \
  https://zenodo.org/record/23326397/files/samtools_1.20.0.img