#!/usr/bin/env python3
"""Summarise the pipeline output in results/tables/ and draw the main figure.

Run from the repository root after the notebook's results have been placed in results/tables/:
    python scripts/summarise_results.py

Writes figures/resistance_mutations.png and results/summary.md.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "results" / "tables"
GENE_ORDER = ["rpoB", "katG", "fabG1", "inhA", "gyrA", "embB", "rpsL", "pncA"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3de", "#fcfcfb"
OK, MISS = "#2a78d6", "#eb6834"


def short(s):
    return s.split(".subj.")[1].split(".lab.")[0] if ".subj." in s else s


def main():
    samples = pd.read_csv(T / "samples.csv")
    qc = pd.read_csv(T / "qc.csv")
    pred = pd.read_csv(T / "predictions.csv").fillna("")
    perf = pd.read_csv(T / "performance.csv")
    v = pd.read_csv(T / "variants_annotated.csv")
    n = samples["sample"].nunique()

    # Variants found in every isolate are differences from the H37Rv reference strain, not resistance.
    v = v.drop_duplicates(["sample", "label"])
    per_label = v.groupby("label")["sample"].nunique()
    shared = sorted(per_label[per_label == n].index)
    v_inf = v[(v.kind != "synonymous") & (~v.label.isin(shared))]

    rule_muts = sorted(set(";".join(pred.RIF_mutations.tolist() + pred.INH_mutations.tolist()).split(";")) - {""})
    other = [l for l in v_inf.label.unique() if l not in rule_muts]
    gene_rank = {g: i for i, g in enumerate(GENE_ORDER)}
    rank = lambda l: (gene_rank.get(l.split("_")[0], 99), l)
    cols = sorted(rule_muts, key=rank) + sorted(other, key=rank)

    grp_order = {"MDR": 0, "INH-mono": 1, "susceptible": 2}
    order = pred.assign(o=pred.group.map(grp_order)).sort_values(["o", "sample"]).reset_index(drop=True)
    mat = np.full((len(order), len(cols)), np.nan)
    for i, s in enumerate(order["sample"]):
        vs = v_inf[v_inf["sample"] == s]
        for j, c in enumerate(cols):
            hit = vs[vs.label == c]
            if len(hit):
                mat[i, j] = hit.alt_fraction.fillna(1).max()

    plt.rcParams.update({"font.size": 8.5, "figure.facecolor": SURF, "savefig.facecolor": SURF, "axes.facecolor": SURF})
    fig, ax = plt.subplots(figsize=(12.5, 5.6))
    for i in range(len(order)):
        for j in range(len(cols)):
            if not np.isnan(mat[i, j]):
                ax.add_patch(plt.Rectangle((j - 0.42, i - 0.38), 0.84, 0.76, color=INK2 if j >= len(rule_muts) else INK,
                                           alpha=0.35 + 0.65 * mat[i, j], lw=0))
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xticks(range(len(cols)), cols, rotation=60, ha="right")
    ax.set_yticks(range(len(order)), [short(s) for s in order["sample"]])
    ax.axvline(len(rule_muts) - 0.5, color=INK2, lw=1)
    for x in np.arange(-0.5, len(cols), 1):
        ax.axvline(x, color=GRID, lw=0.4, zorder=0)
    prev = None
    for i, g in enumerate(order.group):
        if prev is not None and g != prev:
            ax.axhline(i - 0.5, color=INK2, lw=0.8)
        prev = g
    # lab vs predicted, to the right
    for i, r in order.iterrows():
        for k, d in enumerate(("RIF", "INH")):
            lab, pr = r[f"{d}_lab"], r[f"{d}_pred"]
            ax.text(len(cols) + 0.3 + 1.6 * k, i, f"{lab} / {pr}", va="center", fontsize=8.5,
                    color=OK if lab == pr else MISS, fontweight="bold" if lab != pr else "normal", clip_on=False)
    for k, d in enumerate(("RIF", "INH")):
        ax.text(len(cols) + 0.3 + 1.6 * k, -0.9, f"{d}\nlab / pred", fontsize=8, color=INK2, va="bottom", clip_on=False)
    for i, g in order.groupby("group", sort=False):
        ax.text(-0.115, g.index.to_numpy().mean(), i, rotation=90, va="center", ha="center", color=INK2, fontsize=8.5,
                transform=ax.get_yaxis_transform(), clip_on=False)
    ax.text(len(rule_muts) / 2 - 0.5, -1.0, "used by the prediction rules", ha="center", fontsize=8, color=INK)
    ax.text(len(rule_muts) + len(other) / 2 - 0.5, -1.0, "other coding or promoter changes in the panel genes", ha="center", fontsize=8, color=INK2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(length=0)
    fig.suptitle("Drug-resistance mutations called from raw reads in 12 M. tuberculosis isolates", x=0.01, y=0.985, ha="left", fontsize=11, fontweight="bold")
    fig.text(0.01, 0.915, f"Shade = fraction of reads carrying the mutation. Hidden: synonymous changes and {len(shared)} changes present "
             f"in all {n} isolates (differences from the H37Rv reference). Orange = prediction disagrees with the lab.",
             fontsize=8, color=INK2)
    fig.subplots_adjust(left=0.15, right=0.88, top=0.80, bottom=0.27)
    out = ROOT / "figures" / "resistance_mutations.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=200)

    # ---- summary.md ----
    q = qc.set_index("sample")
    L = ["# Results summary", "", "Generated by `scripts/summarise_results.py` from `results/tables/`.", "",
         "## Sequencing quality", "",
         f"- Median depth across the genome: {qc.mean_depth.median():.0f}x (range {qc.mean_depth.min():.0f}-{qc.mean_depth.max():.0f}x).",
         f"- Genome covered: {qc.genome_covered_pct.min():.1f}-{qc.genome_covered_pct.max():.1f}% "
         f"(all but one isolate at least {qc.genome_covered_pct.nlargest(n - 1).min():.2f}%).",
         f"- Minimum depth in rpoB, katG and fabG1 across isolates: {qc[['depth_rpoB', 'depth_katG', 'depth_fabG1']].min().min():.0f}x.", "",
         "## Agreement with laboratory drug-susceptibility testing", "",
         "| Drug | Resistant in lab, called R | Resistant in lab, called S | Susceptible in lab, called S | Susceptible in lab, called R | Sensitivity (95% CI) | Specificity (95% CI) |",
         "|---|---|---|---|---|---|---|"]
    for _, r in perf.iterrows():
        L.append(f"| {r.drug} | {r.TP} | {r.FN} | {r.TN} | {r.FP} | {r.sensitivity:.0%} ({r.sens_lo:.0%}-{r.sens_hi:.0%}) "
                 f"| {r.specificity:.0%} ({r.spec_lo:.0%}-{r.spec_hi:.0%}) |")
    L += ["", "## Mutations used by the rules", ""]
    for m in rule_muts:
        L.append(f"- {m}: {int((v.label == m).groupby(v['sample']).any().sum())} isolate(s)")
    L += ["", f"## Changes present in all {n} isolates (reference differences, hidden from the figure)", "", "- " + ", ".join(shared), "",
          "## Discordant isolates", ""]
    for d in ("RIF", "INH"):
        for _, r in pred[pred[f"{d}_lab"] != pred[f"{d}_pred"]].iterrows():
            nonsyn = v_inf[v_inf["sample"] == r["sample"]].label.tolist()
            L.append(f"- {d}, {short(r['sample'])} ({r.group}): lab {r[f'{d}_lab']}, predicted {r[f'{d}_pred']}. "
                     f"Coding or promoter changes in panel genes (excluding reference differences): {', '.join(nonsyn) if nonsyn else 'none'}.")
    (ROOT / "results" / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
