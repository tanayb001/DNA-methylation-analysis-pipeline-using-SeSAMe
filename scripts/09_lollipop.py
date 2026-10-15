"""STEP 9 - Lollipop plots for every pathway table (FDR < cutoff).

If a table has no term below the FDR cutoff, the top terms by NOMINAL p-value are
drawn instead and the figure is clearly labelled 'EXPLORATORY' - never as significant.
"""
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
res = Path(cfg["paths"]["results_dir"])
tab, out = res / "tables", res / "figures" / "pathways"
out.mkdir(parents=True, exist_ok=True)
fdr, top_n = cfg["pathway"]["fdr_cutoff"], cfg["pathway"]["top_n_plot"]

for f in sorted(tab.glob("pathway_*.csv")):
    df = pd.read_csv(f)
    if df.empty:
        continue
    gsea = "NES" in df.columns
    sig = df[df["p.adjust"] < fdr]
    exploratory = sig.empty
    if exploratory:
        sig = df.nsmallest(top_n, "pvalue")
    else:
        sig = sig.nsmallest(top_n, "p.adjust")
    sig = sig.copy()

    if gsea:
        x, xlabel = sig["NES"].values, "Normalised enrichment score (NES)"
        size = sig["setSize"].values
    else:
        col = "pvalue" if exploratory else "p.adjust"
        x = -np.log10(sig[col].values)
        xlabel = f"-log10({'nominal p' if exploratory else 'FDR'})"
        size = sig["Count"].values
    colour = -np.log10(sig["pvalue" if exploratory else "p.adjust"].values)
    order = np.argsort(x)
    sig, x, size, colour = sig.iloc[order], x[order], size[order], colour[order]
    labels = ["\n".join(textwrap.wrap(str(t).replace("_", " "), 55)) for t in sig["Description"]]

    fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.42 * len(sig) + 1.5)))
    ypos = np.arange(len(sig))
    ax.hlines(ypos, 0, x, color="grey", linewidth=1.2)
    sc = ax.scatter(x, ypos, s=40 + 12 * (size / max(size.max(), 1)) * 10, c=colour,
                    cmap="viridis_r", edgecolor="black", zorder=3)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(xlabel)
    title = f.stem.replace("pathway_", "")
    ax.set_title(title + ("\nEXPLORATORY: no term passed FDR < %s (top by nominal p)" % fdr if exploratory
                          else f"\nFDR < {fdr}"), fontsize=10,
                 color="firebrick" if exploratory else "black")
    cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label("-log10(" + ("p" if exploratory else "FDR") + ")")
    fig.tight_layout()
    fig.savefig(out / f"lollipop_{title}.pdf"); fig.savefig(out / f"lollipop_{title}.png", dpi=300)
    plt.close(fig)
    print(f"{title}: {'EXPLORATORY' if exploratory else 'FDR-significant'} ({len(sig)} terms)")
