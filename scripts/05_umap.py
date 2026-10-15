"""STEP 5 - UMAP of samples using the most variable probes (exported by step 4)."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import umap
import yaml

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
res = Path(cfg["paths"]["results_dir"])
tab, fig = res / "tables", res / "figures"

X = pd.read_csv(tab / "top_variable_betas.csv", index_col=0)      # probes x units
ph = pd.read_csv(tab / "pheno_analysis.csv", index_col=0).loc[X.columns]
n = X.shape[1]
if n < 4:
    sys.exit("Fewer than 4 samples: UMAP is not meaningful, skipping.")

nn = int(max(2, min(cfg["exploratory"]["umap_n_neighbors"], n - 1)))
reducer = umap.UMAP(n_neighbors=nn, min_dist=cfg["exploratory"]["umap_min_dist"],
                    metric="correlation", random_state=cfg["seed"],
                    init="random" if n < 15 else "spectral")
emb = reducer.fit_transform(X.T.values)
out = pd.DataFrame(emb, index=X.columns, columns=["UMAP1", "UMAP2"])
out["group"] = ph["group"].values
out.to_csv(tab / "umap_coordinates.csv")

palette = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
fig_, ax = plt.subplots(figsize=(6, 5))
for i, (g, sub) in enumerate(out.groupby("group")):
    ax.scatter(sub.UMAP1, sub.UMAP2, s=70, color=palette[i % len(palette)], label=g)
for name, row in out.iterrows():
    ax.annotate(name, (row.UMAP1, row.UMAP2), fontsize=8, xytext=(4, 4), textcoords="offset points")
ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
ax.set_title(f"UMAP - top {X.shape[0]} variable probes (n_neighbors={nn}, n={n})")
ax.legend(title="Group", frameon=False)
fig_.tight_layout()
fig_.savefig(fig / "06_UMAP.pdf"); fig_.savefig(fig / "06_UMAP.png", dpi=300)
print("UMAP done (note: with very few samples treat UMAP as descriptive only).")
