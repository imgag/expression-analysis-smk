import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

metadata = pd.read_csv(snakemake.input["metadata"], sep="\t", index_col=0)
contrasts = pd.read_csv(snakemake.input["contrasts"], sep="\t", index_col=0)
dat = pd.read_csv(snakemake.input["contrast"], sep="\t")

dat["logFDR"] = -np.log10(dat["FDR"])
dat["significant"] = pd.Categorical(dat["FDR"] < snakemake.config["de_filter"]["fdr"])
contrast = snakemake.wildcards["contrast"]
p = dat.plot.scatter(
    x="logFC", y="logFDR", c=dat["significant"].map({True: "red", False: "grey"})
)
p.set(xlabel="log2 FC", ylabel="-log10 FDR", title=contrasts.loc[contrast, "title"])
plt.savefig(snakemake.output["plot"], dpi=300, bbox_inches="tight")
