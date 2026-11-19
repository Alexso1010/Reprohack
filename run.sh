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
# --- Suivi de l'exécution ---------------------------------------------------
mkdir -p logs
LOG="logs/run_$(date +%Y%m%d_%H%M%S).log"
HEARTBEAT="${HEARTBEAT:-60}"     # intervalle en secondes entre deux points

heartbeat() {
    while true; do
        sleep "${HEARTBEAT}"
        progress=$(grep -o '[0-9]* of [0-9]* steps ([0-9]*%) done' "$LOG" | tail -n 1)
        echo "[$(date +%H:%M:%S)] ${progress:-démarrage...}" \
             "| charge CPU: $(cut -d' ' -f1 /proc/loadavg)" \
             "| RAM dispo: $(free -h | awk '/Mem:/{print $7}')" \
             "| disque libre: $(df -h . | awk 'NR==2{print $4}')"
    done
}
heartbeat &
HB_PID=$!
trap 'kill "$HB_PID" 2>/dev/null' EXIT     # arrête le suivi quand le script se termine

# Toute la sortie de Snakemake est affichée ET enregistrée dans $LOG
snakemake --use-apptainer --cores "${CORES}" --printshellcmds "$@" 2>&1 | tee -a "$LOG"