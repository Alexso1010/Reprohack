#!/bin/bash
# Relance toute l'analyse du projet 1 (persisters de S. aureus).
#
# Usage :
#   bash run.sh            # exécution complète, 4 cœurs
#   CORES=8 bash run.sh    # choisir le nombre de cœurs
#   bash run.sh -n         # "dry-run" : vérifie le workflow sans rien exécuter
#
# Prérequis : apptainer et snakemake installés (voir README.md).

set -euo pipefail #arrêt du scrip si erreur

# Se placer à la racine du dépôt, quel que soit l'endroit d'où on lance le script
cd "$(dirname "$0")"

CORES="${CORES:-4}"

# 1. Vérifier que les outils nécessaires sont présents
for cmd in apptainer snakemake wget; do
    if ! command -v "$cmd" > /dev/null 2>&1; then   #command -v vérifie qu'on a bien Snakemake et Apptainer
        echo "Erreur : '$cmd' est introuvable. Voir les prérequis dans le README.md." >&2
        exit 1
    fi
done

# 2. Préparer les conteneurs (téléchargement Damona + construction des recettes .def)
bash build_containers.sh

# 3. Lancer le workflow ; les arguments éventuels (ex. -n) sont transmis à Snakemake
snakemake --use-apptainer --cores "${CORES}" "$@"