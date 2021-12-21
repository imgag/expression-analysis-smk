import matplotlib.cm as cm
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.patches import Patch

matplotlib.use("Agg")
import matplotlib.pyplot as plt

metadata = pd.read_csv(snakemake.input["metadata"], sep="\t", index_col="sample")
contrast = pd.read_csv(snakemake.input["contrast"], sep="\t", index_col=0)
expr = pd.read_csv(snakemake.input["filtered_tpm"], sep="\t", index_col=0).drop(
    columns=["gene_length"]
)

genes = contrast.head(snakemake.config["heatmap_n_genes"]).index
expr_plot = (
    expr.loc[genes, :]
    .merge(contrast[["symbol"]], left_index=True, right_index=True, how="left")
    .set_index("symbol")
)

group_colors_df = pd.read_csv(snakemake.input["group_colors"], sep="\t", index_col=0)
group_colors_dict = dict(group_colors_df["color"])
row_colors = list(metadata["group"].map(group_colors_dict))

p = sns.clustermap(
    expr_plot.apply(lambda x: np.log2(x + 1)),
    method="average",
    col_colors=row_colors,
    cmap=snakemake.config["colormaps"]["heatmap"],
    z_score=0,
    yticklabels=True,
)
fig = plt.gcf()
fig.set_size_inches(8, 20)

handles = [Patch(facecolor=c) for c in group_colors_dict.values()]
plt.legend(
    handles,
    group_colors_dict.keys(),
    bbox_to_anchor=(1, 1),
    bbox_transform=plt.gcf().transFigure,
    loc="upper right",
)

plt.savefig(snakemake.output["heatmap"], dpi=300, bbox_inches="tight")
