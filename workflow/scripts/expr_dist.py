import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# d_full_cpm = pd.read_csv(snakemake.input["full_cpm"], sep="\t", index_col=[0, 1])
# d_full_cpm_log2 = d_full_cpm.apply(lambda x: np.log2(x + 1))

# ax = d_full_cpm_log2.boxplot(rot=90, grid=False)
# ax.set(title="Distribution of unfiltered expression values", ylabel="log2(CPM+1)")
# plt.savefig(snakemake.output["boxplot_full"], dpi=300, bbox_inches="tight")

# ax = d_full_cpm_log2.plot.kde()
# ax.set(title="Distribution of unfiltered expression values", xlabel="log2(CPM+1)")
# plt.savefig(snakemake.output["density_full"], dpi=300, bbox_inches="tight")

d_cpm = pd.read_csv(snakemake.input["filtered_cpm"], sep="\t", index_col=[0, 1])
d_cpm_log2 = d_cpm.apply(lambda x: np.log2(x + 1))

ax = d_cpm_log2.boxplot(rot=90, grid=False)
ax.set(title="Distribution of filtered expression values", ylabel="log2(CPM+1)")
plt.savefig(snakemake.output["boxplot"], dpi=300, bbox_inches="tight")

ax = d_cpm_log2.plot.kde()
ax.set(title="Distribution of filtered expression values", xlabel="log2(CPM+1)")
plt.savefig(snakemake.output["density"], dpi=300, bbox_inches="tight")
