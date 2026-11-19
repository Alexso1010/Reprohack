#!/usr/bin/env python3
"""Compare deux tables de résultats DESeq2 : une référence (celle des auteurs de
l'article) et notre sortie.

Exemple :
    python3 compare_results.py reference.csv resultats/deseq2_results.csv -o comparaison/

Colonnes attendues
    référence : Id, Name, baseMean, log2FoldChange, pvalue, padj (+ autres, ignorées)
    notre sortie : baseMean, log2FoldChange, pvalue, padj, locus, name (+ autres, ignorées)

Les deux fichiers sont mis en correspondance via l'identifiant du gène :
    colonne `Id` de la référence  <->  colonne `locus` de notre sortie
(modifiable avec --ref-id et --ours-id).

Le séparateur (virgule, tabulation, point-virgule) est détecté automatiquement.

Sorties (dans le dossier donné par -o) :
    summary.txt            rapport lisible (aussi affiché à l'écran)
    metrics.json           métriques sous forme exploitable par un programme
    merged.csv             table fusionnée, une ligne par gène (union des deux fichiers)
    discordant_genes.csv   gènes dont la significativité ou le sens diffère
    top_differences.csv    gènes avec les plus grands écarts de log2FoldChange
    comparison_plots.png   figures (si matplotlib est installé)

Dépendances : pandas, numpy (matplotlib facultatif, pour les figures).
"""
/Users/annagarnon/Documents/AGROPARISTECH/3AIODAA/Reproductibilite_computationnelle/Reprohack/deseq2_results.csv


import argparse
import json
import sys
from pathlib import Path


import numpy as np
import pandas as pd

# Colonnes utilisées dans chaque fichier
REF_COLS = {"name": "Name", "baseMean": "baseMean", "lfc": "log2FoldChange",
            "pvalue": "pvalue", "padj": "padj"}
OURS_COLS = {"name": "name", "baseMean": "baseMean", "lfc": "log2FoldChange",
             "pvalue": "pvalue", "padj": "padj"}


# ----------------------------------------------------------------------------
# Lecture des fichiers
# ----------------------------------------------------------------------------
def detect_sep(path):
    """Devine le séparateur d'après la ligne d'en-tête."""
    with open(path, encoding="utf-8-sig") as handle:
        first = handle.readline()
    counts = {sep: first.count(sep) for sep in (",", "\t", ";")}
    return max(counts, key=counts.get)


def normalize_id(value):
    """Uniformise les identifiants : espaces retirés, préfixe 'gene-' (GFF NCBI) retiré."""
    if pd.isna(value):
        return value
    value = str(value).strip().strip('"')
    if value.startswith("gene-"):
        value = value[len("gene-"):]
    return value


def read_table(path, id_col, cols, label, sep=None):
    """Lit un fichier, vérifie les colonnes, renvoie un DataFrame indexé par identifiant normalisé."""
    sep = sep or detect_sep(path)
    df = pd.read_csv(path, sep=sep, encoding="utf-8-sig")
    needed = [id_col] + list(cols.values())
    missing = [c for c in needed if c not in df.columns]
    if missing:
        sys.exit(f"Erreur : colonnes absentes du fichier {label} ({path}) : {missing}\n"
                 f"Colonnes trouvées : {list(df.columns)}")
    out = pd.DataFrame({
        "id": df[id_col].map(normalize_id),
        "name": df[cols["name"]],
        "baseMean": pd.to_numeric(df[cols["baseMean"]], errors="coerce"),
        "log2FC": pd.to_numeric(df[cols["lfc"]], errors="coerce"),
        "pvalue": pd.to_numeric(df[cols["pvalue"]], errors="coerce"),
        "padj": pd.to_numeric(df[cols["padj"]], errors="coerce"),
    })
    return out


# ----------------------------------------------------------------------------
# Petites fonctions statistiques (sans scipy)
# ----------------------------------------------------------------------------
def pearson(x, y):
    mask = x.notna() & y.notna()
    if mask.sum() < 3:
        return float("nan")
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def spearman(x, y):
    mask = x.notna() & y.notna()
    if mask.sum() < 3:
        return float("nan")
    return pearson(x[mask].rank(), y[mask].rank())


def linear_fit(x, y):
    """Pente et ordonnée à l'origine de y = a*x + b (x = référence, y = nous)."""
    mask = x.notna() & y.notna()
    if mask.sum() < 3:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(x[mask], y[mask], 1)
    return float(slope), float(intercept)


def cohen_kappa(both, ref_only, ours_only, neither):
    n = both + ref_only + ours_only + neither
    if n == 0:
        return float("nan")
    observed = (both + neither) / n
    expected = ((both + ref_only) * (both + ours_only)
                + (ours_only + neither) * (ref_only + neither)) / n ** 2
    return float((observed - expected) / (1 - expected)) if expected < 1 else float("nan")


def ratio(a, b):
    return float(a / b) if b else float("nan")


def significance(df, suffix, alpha, lfc_min):
    """Booléens : significatif / up / down. padj ou log2FC manquant => non significatif."""
    padj = df[f"padj_{suffix}"]
    lfc = df[f"log2FC_{suffix}"]
    sig = (padj < alpha) & (lfc.abs() >= lfc_min)
    return sig.fillna(False), (sig & (lfc > 0)).fillna(False), (sig & (lfc < 0)).fillna(False)


# ----------------------------------------------------------------------------
# Métriques
# ----------------------------------------------------------------------------
def file_qc(df, label, alpha, lfc_min):
    """Contrôles sur un fichier pris isolément."""
    sig = (df["padj"] < alpha) & (df["log2FC"].abs() >= lfc_min)
    return {
        "fichier": label,
        "lignes": int(len(df)),
        "identifiants_uniques": int(df["id"].nunique()),
        "identifiants_manquants": int(df["id"].isna().sum()),
        "log2FC_NA": int(df["log2FC"].isna().sum()),
        "padj_NA": int(df["padj"].isna().sum()),
        "significatifs": int(sig.sum()),
        "up": int((sig & (df["log2FC"] > 0)).sum()),
        "down": int((sig & (df["log2FC"] < 0)).sum()),
    }


def compare(merged, alpha, lfc_min):
    """Métriques de comparaison sur les gènes présents dans les deux fichiers."""
    common = merged[merged["_merge"] == "both"].copy()
    m = {"genes_communs": int(len(common))}
    if common.empty:
        return m, common

    # Corrélations et écarts
    m["pearson_log2FC"] = pearson(common["log2FC_ref"], common["log2FC_ours"])
    m["spearman_log2FC"] = spearman(common["log2FC_ref"], common["log2FC_ours"])
    slope, intercept = linear_fit(common["log2FC_ref"], common["log2FC_ours"])
    m["pente_log2FC_(nous=a*ref+b)"] = slope
    m["ordonnee_log2FC"] = intercept
    delta = common["log2FC_ours"] - common["log2FC_ref"]
    m["ecart_moyen_log2FC_(nous-ref)"] = float(delta.mean())
    m["ecart_absolu_moyen_log2FC"] = float(delta.abs().mean())
    m["RMSE_log2FC"] = float(np.sqrt((delta ** 2).mean()))
    m["gènes_avec_ecart_log2FC>0.5"] = int((delta.abs() > 0.5).sum())
    m["gènes_avec_ecart_log2FC>1"] = int((delta.abs() > 1).sum())

    log_ref = np.log2(common["baseMean_ref"] + 1)
    log_ours = np.log2(common["baseMean_ours"] + 1)
    m["pearson_log2(baseMean+1)"] = pearson(log_ref, log_ours)
    m["spearman_baseMean"] = spearman(common["baseMean_ref"], common["baseMean_ours"])
    m["ratio_median_baseMean_(nous/ref)"] = float(
        (common["baseMean_ours"] / common["baseMean_ref"]).replace([np.inf, -np.inf], np.nan).median())

    # Spearman sur -log10(padj) : mesure l'accord sur la force du signal
    nlp_ref = -np.log10(common["padj_ref"].clip(lower=1e-300))
    nlp_ours = -np.log10(common["padj_ours"].clip(lower=1e-300))
    m["spearman_-log10(padj)"] = spearman(nlp_ref, nlp_ours)

    # Accord sur le sens (signe du log2FC), tous gènes communs
    m["accord_signe_log2FC_%"] = float(
        100 * (np.sign(common["log2FC_ref"]) == np.sign(common["log2FC_ours"])).mean())

    # Accord sur les gènes différentiellement exprimés (DEG)
    sig_r, up_r, down_r = significance(common, "ref", alpha, lfc_min)
    sig_o, up_o, down_o = significance(common, "ours", alpha, lfc_min)
    both = int((sig_r & sig_o).sum())
    ref_only = int((sig_r & ~sig_o).sum())
    ours_only = int((~sig_r & sig_o).sum())
    neither = int((~sig_r & ~sig_o).sum())
    m["DEG_ref_(parmi_communs)"] = int(sig_r.sum())
    m["DEG_nous_(parmi_communs)"] = int(sig_o.sum())
    m["DEG_dans_les_deux"] = both
    m["DEG_ref_seulement"] = ref_only
    m["DEG_nous_seulement"] = ours_only
    m["DEG_ni_l_un_ni_l_autre"] = neither
    m["rappel_(DEG_ref_retrouvés)"] = ratio(both, both + ref_only)
    m["precision_(nos_DEG_dans_ref)"] = ratio(both, both + ours_only)
    m["jaccard_DEG"] = ratio(both, both + ref_only + ours_only)
    m["kappa_de_Cohen_DEG"] = cohen_kappa(both, ref_only, ours_only, neither)

    both_mask = sig_r & sig_o
    same_dir = ((up_r & up_o) | (down_r & down_o))[both_mask]
    m["DEG_communs_meme_sens"] = int(same_dir.sum())
    m["DEG_communs_sens_oppose"] = int(both_mask.sum() - same_dir.sum())
    m["up_ref_et_nous"] = int((up_r & up_o).sum())
    m["down_ref_et_nous"] = int((down_r & down_o).sum())

    # Catégorie par gène (utile pour les figures et la liste des discordances)
    category = np.select(
        [sig_r & sig_o, sig_r & ~sig_o, ~sig_r & sig_o],
        ["DEG dans les deux", "DEG référence seulement", "DEG nous seulement"],
        default="non significatif")
    common["categorie"] = category
    common["delta_log2FC"] = delta
    common["sens_oppose"] = both_mask & ~((up_r & up_o) | (down_r & down_o))
    return m, common


# ----------------------------------------------------------------------------
# Mise en forme du rapport
# ----------------------------------------------------------------------------
def fmt(value):
    if isinstance(value, float):
        return "NA" if np.isnan(value) else f"{value:.4g}"
    return str(value)


def build_report(qc_ref, qc_ours, overlap, metrics, title, alpha, lfc_min):
    lines = [f"===== {title} =====", ""]
    lines.append(f"Seuils DEG utilisés : padj < {alpha}, |log2FC| >= {lfc_min}")
    lines.append("")
    lines.append("-- Contrôle de chaque fichier --")
    for qc in (qc_ref, qc_ours):
        lines.append(f"[{qc['fichier']}]")
        for key, val in qc.items():
            if key != "fichier":
                lines.append(f"  {key:28s} {fmt(val)}")
    lines.append("")
    lines.append("-- Recouvrement des gènes --")
    for key, val in overlap.items():
        lines.append(f"  {key:28s} {fmt(val)}")
    lines.append("")
    lines.append("-- Comparaison sur les gènes communs --")
    for key, val in metrics.items():
        lines.append(f"  {key:34s} {fmt(val)}")
    lines.append("")
    return lines


def warnings_for(metrics, overlap):
    notes = []
    if overlap["communs"] == 0:
        notes.append("AUCUN gène commun : vérifier --ref-id / --ours-id et le format des identifiants.")
    elif overlap["communs"] < 0.5 * min(overlap["dans_ref"], overlap["dans_nous"]):
        notes.append("Moins de la moitié des gènes sont communs : probable décalage d'identifiants.")
    p = metrics.get("pearson_log2FC", float("nan"))
    if not np.isnan(p) and p < 0:
        notes.append("Corrélation NÉGATIVE des log2FC : le contraste est probablement inversé "
                     "(persisters/contrôle vs contrôle/persisters).")
    if metrics.get("DEG_communs_sens_oppose", 0) > 0:
        notes.append(f"{metrics['DEG_communs_sens_oppose']} gènes significatifs des deux côtés "
                     "mais dans des sens opposés.")
    return notes


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def make_plots(common, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib absent : figures non générées.", file=sys.stderr)
        return False

    colors = {"DEG dans les deux": "tab:red", "DEG référence seulement": "tab:blue",
              "DEG nous seulement": "tab:orange", "non significatif": "lightgrey"}
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    ax = axes[0]
    for cat, col in colors.items():
        sub = common[common["categorie"] == cat]
        ax.scatter(sub["log2FC_ref"], sub["log2FC_ours"], s=8, c=col, label=f"{cat} ({len(sub)})")
    lim = np.nanmax(np.abs(common[["log2FC_ref", "log2FC_ours"]].to_numpy()))
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8)
    ax.set_xlabel("log2FC référence")
    ax.set_ylabel("log2FC nous")
    ax.set_title("Comparaison des log2FoldChange")
    ax.legend(fontsize=7)

    ax = axes[1]
    x = np.log2(common["baseMean_ref"] + 1)
    y = np.log2(common["baseMean_ours"] + 1)
    ax.scatter(x, y, s=8, c="tab:green")
    top = np.nanmax([x.max(), y.max()])
    ax.plot([0, top], [0, top], "k--", lw=0.8)
    ax.set_xlabel("log2(baseMean + 1) référence")
    ax.set_ylabel("log2(baseMean + 1) nous")
    ax.set_title("Niveau d'expression moyen")

    ax = axes[2]
    ax.hist(common["delta_log2FC"].dropna(), bins=60, color="tab:purple")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("log2FC nous - log2FC référence")
    ax.set_ylabel("nombre de gènes")
    ax.set_title("Distribution des écarts")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


# ----------------------------------------------------------------------------
# Programme principal
# ----------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Compare une table DESeq2 de référence à notre table DESeq2.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("reference", help="fichier de référence (résultats des auteurs)")
    p.add_argument("ours", help="notre fichier de résultats")
    p.add_argument("-o", "--outdir", default="comparison", help="dossier de sortie")
    p.add_argument("--ref-id", default="Id", help="colonne identifiant dans la référence")
    p.add_argument("--ours-id", default="locus", help="colonne identifiant dans notre fichier")
    p.add_argument("--alpha", type=float, default=0.05, help="seuil de padj pour les DEG")
    p.add_argument("--lfc", type=float, default=0.0,
                   help="seuil de |log2FC| pour les DEG (0 = pas de seuil, comme dans l'article)")
    p.add_argument("--subset", help="fichier de locus tags (un par ligne) : "
                                    "répète la comparaison sur ce sous-ensemble (ex. gènes de la traduction)")
    p.add_argument("--top", type=int, default=20, help="nombre de gènes dans la liste des plus gros écarts")
    p.add_argument("--sep-ref", help="séparateur de la référence (auto par défaut)")
    p.add_argument("--sep-ours", help="séparateur de notre fichier (auto par défaut)")
    p.add_argument("--no-plots", action="store_true", help="ne pas générer les figures")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ref = read_table(args.reference, args.ref_id, REF_COLS, "référence", args.sep_ref)
    ours = read_table(args.ours, args.ours_id, OURS_COLS, "nous", args.sep_ours)

    notes = []
    for label, df in (("référence", ref), ("nous", ours)):
        dup = int(df["id"].duplicated().sum())
        if dup:
            notes.append(f"{dup} identifiants dupliqués dans {label} : seule la première occurrence est gardée.")
        df.drop_duplicates("id", keep="first", inplace=True)
        df.dropna(subset=["id"], inplace=True)

    qc_ref = file_qc(ref, "référence", args.alpha, args.lfc)
    qc_ours = file_qc(ours, "nous", args.alpha, args.lfc)

    merged = ref.merge(ours, on="id", how="outer", suffixes=("_ref", "_ours"), indicator=True)
    merged["_merge"] = merged["_merge"].astype(str).replace(
        {"left_only": "ref_seulement", "right_only": "nous_seulement"})
    overlap = {
        "dans_ref": int(len(ref)),
        "dans_nous": int(len(ours)),
        "communs": int((merged["_merge"] == "both").sum()),
        "ref_seulement": int((merged["_merge"] == "ref_seulement").sum()),
        "nous_seulement": int((merged["_merge"] == "nous_seulement").sum()),
    }

    metrics, common = compare(merged, args.alpha, args.lfc)
    notes += warnings_for(metrics, overlap)

    lines = build_report(qc_ref, qc_ours, overlap, metrics,
                         "COMPARAISON GLOBALE", args.alpha, args.lfc)
    all_metrics = {"qc_reference": qc_ref, "qc_nous": qc_ours, "recouvrement": overlap,
                   "comparaison": metrics}

    # Sous-ensemble facultatif (ex. gènes de la traduction)
    if args.subset:
        subset_ids = {normalize_id(x) for x in Path(args.subset).read_text().split()}
        sub_merged = merged[merged["id"].isin(subset_ids)]
        sub_overlap = {
            "dans_le_sous_ensemble": len(subset_ids),
            "dans_ref": int(sub_merged["log2FC_ref"].notna().sum()),
            "dans_nous": int(sub_merged["log2FC_ours"].notna().sum()),
            "communs": int((sub_merged["_merge"] == "both").sum()),
        }
        sub_metrics, _ = compare(sub_merged, args.alpha, args.lfc)
        lines += ["", f"===== SOUS-ENSEMBLE : {args.subset} =====", "", "-- Recouvrement --"]
        lines += [f"  {k:28s} {fmt(v)}" for k, v in sub_overlap.items()]
        lines += ["", "-- Comparaison sur les gènes communs du sous-ensemble --"]
        lines += [f"  {k:34s} {fmt(v)}" for k, v in sub_metrics.items()]
        all_metrics["sous_ensemble"] = {"recouvrement": sub_overlap, "comparaison": sub_metrics}

    if notes:
        lines += ["", "-- À SURVEILLER --"] + [f"  ! {n}" for n in notes]

    # Fichiers de sortie
    merged.to_csv(outdir / "merged.csv", index=False)
    if not common.empty:
        discordant = common[(common["categorie"].isin(["DEG référence seulement", "DEG nous seulement"]))
                            | common["sens_oppose"]]
        discordant = discordant.sort_values("delta_log2FC", key=lambda s: s.abs(), ascending=False)
        discordant.to_csv(outdir / "discordant_genes.csv", index=False)
        top = common.sort_values("delta_log2FC", key=lambda s: s.abs(), ascending=False).head(args.top)
        top.to_csv(outdir / "top_differences.csv", index=False)
        if not args.no_plots:
            if make_plots(common, outdir / "comparison_plots.png"):
                lines.append(f"\nFigures : {outdir / 'comparison_plots.png'}")

    report = "\n".join(lines)
    (outdir / "summary.txt").write_text(report + "\n", encoding="utf-8")
    with open(outdir / "metrics.json", "w", encoding="utf-8") as handle:
        json.dump(all_metrics, handle, indent=2, ensure_ascii=False, default=float)
    print(report)
    print(f"\nRésultats écrits dans : {outdir}/")


if __name__ == "__main__":
    main()
