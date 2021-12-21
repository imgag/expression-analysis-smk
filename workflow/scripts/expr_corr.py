import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

d = pd.read_csv(snakemake.input["filtered_raw"], sep="\t", index_col=[0, 1])

corr = d.corr(method="spearman").abs()
corr.rename_axis("sample").reset_index().to_csv(
    snakemake.output["table"], sep="\t", index=False
)

linked = linkage(squareform(1 - corr), method="average")

f, ax = plt.subplots()
dendrogram(
    linked,
    orientation="left",
    labels=corr.index,
    distance_sort="descending",
    show_leaf_counts=True,
    ax=ax,
)
ax.set(
    title="Hierarchical clustering on correlation-based distance",
    xlabel="1 - rank correlation coefficient",
)
plt.savefig(snakemake.output["dendrogram"], dpi=300, bbox_inches="tight")
