import matplotlib.cm as cm
import seaborn as sns
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

metadata = pd.read_csv(snakemake.input["metadata"], sep="\t", index_col="sample")
corr = pd.read_csv(snakemake.input["corr_tbl"], sep="\t", index_col=0)

group_colors_df = pd.read_csv(snakemake.input["group_colors"], sep="\t", index_col=0)
group_colors_dict = dict(group_colors_df["color"])
row_colors = list(metadata["group"].map(group_colors_dict))
sns.clustermap(
    corr,
    method="average",
    row_colors=row_colors,
    col_colors=row_colors,
    vmin=0,
    vmax=1,
    cmap=snakemake.config["colormaps"]["clustermap"],
    xticklabels=True,
    yticklabels=True,
)

handles = [Patch(facecolor=c) for c in group_colors_dict.values()]
plt.legend(
    handles,
    group_colors_dict.keys(),
    bbox_to_anchor=(1, 1),
    bbox_transform=plt.gcf().transFigure,
    loc="upper right",
)

plt.savefig(snakemake.output["clustermap"], dpi=300, bbox_inches="tight")
