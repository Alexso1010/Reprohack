#!/usr/bin/env python3
"""Liste les gènes (locus tags) d'une ou plusieurs voies KEGG, via l'API REST.

Usage : python3 kegg_translation_genes.py sao03010 sao00970 > translation_genes.txt

Sortie : un locus tag par ligne (ex. SAOUHSC_00481), dédoublonné et trié.
"""
import sys
import urllib.request


def parse_link(text):
    """Extrait les locus tags d'une réponse KEGG 'link' (lignes 'path:sao03010<TAB>sao:SAOUHSC_00481')."""
    genes = set()
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) == 2 and parts[1].startswith("sao:"):
            genes.add(parts[1][len("sao:"):])
    return genes


def fetch_pathway_genes(pathway):
    url = f"https://rest.kegg.jp/link/sao/{pathway}"
    with urllib.request.urlopen(url, timeout=60) as response:
        return parse_link(response.read().decode("utf-8"))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    all_genes = set()
    for pathway in sys.argv[1:]:
        found = fetch_pathway_genes(pathway)
        if not found:
            sys.exit(f"Erreur : aucun gène renvoyé par KEGG pour {pathway}")
        all_genes |= found
    print("\n".join(sorted(all_genes)))
