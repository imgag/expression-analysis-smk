import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

n_components = 3

metadata = pd.read_csv(snakemake.input["metadata"], sep="\t", index_col=0)
dat = pd.read_csv(snakemake.input["filtered_cpm"], sep="\t")

X = dat.iloc[:, 2:].T
pca = PCA(n_components=n_components)
pca.fit(X)
X_pca = pca.transform(X)
X_pca_df = pd.DataFrame(X_pca, index=X.index).merge(
    metadata, left_index=True, right_index=True, how="left"
)
X_pca_df["sample"] = X_pca_df.index

group_colors_df = pd.read_csv(snakemake.input["group_colors"], sep="\t", index_col=0)
group_colors_dict = dict(group_colors_df["color"])
p = sns.scatterplot(
    x=0,
    y=1,
    data=X_pca_df,
    hue=snakemake.config["pca"]["hue"],
    palette=group_colors_dict,
    style=snakemake.config["pca"]["style"],
    alpha=0.75,
)
p.set(
    xlabel="PC1: {:2.2%}".format(pca.explained_variance_ratio_[0]),
    ylabel="PC2: {:2.2%}".format(pca.explained_variance_ratio_[1]),
)
plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
plt.axis("equal")
if snakemake.config["pca"]["label"]:
    for idx, row in X_pca_df.iterrows():
        plt.text(row[0], row[1], row[snakemake.config["pca"]["label"]], fontsize=6, color=group_colors_dict[row[snakemake.config["pca"]["hue"]]], horizontalalignment="center")

plt.savefig(snakemake.output["pca_plt"], dpi=300, bbox_inches="tight")
X_pca_df.iloc[:, 0:n_components].to_csv(snakemake.output["pca_tbl"], sep="\t")
