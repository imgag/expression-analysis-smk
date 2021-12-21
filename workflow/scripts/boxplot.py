import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


metadata = pd.read_csv(snakemake.input["metadata"], sep="\t")
dat = pd.read_csv(snakemake.input["filtered_tpm_annot"], sep="\t")
gene_symbol = snakemake.wildcards["gene"]

dat = dat[dat["symbol"] == gene_symbol]

assert len(dat.index) == 1

dat_long = dat.melt(
    id_vars=["gene_id", "symbol", "description", "biotype", "gene_length"],
    var_name="sample",
    value_name="tpm",
)
dat_long_meta = dat_long.merge(
    metadata[["sample", "group"]], on="sample", validate="m:1"
)
dat_long_meta["logTPM"] = dat_long_meta["tpm"].apply(lambda x: np.log2(x + 1))

f, ax = plt.subplots()
dat_long_meta[["sample", "logTPM", "group", "gene_id", "symbol"]].groupby(
    ["gene_id"]
).boxplot(column="logTPM", ax=ax, by="group")
ax.set(title=gene_symbol, xlabel="expression", ylabel="log2(TPM+1)")
plt.suptitle("")
plt.savefig(snakemake.output["boxplot_log"], dpi=300, bbox_inches="tight")

f, ax = plt.subplots()
dat_long_meta[["sample", "tpm", "group", "gene_id", "symbol"]].groupby(
    ["gene_id"]
).boxplot(column="tpm", ax=ax, by="group")
ax.set(title=gene_symbol, xlabel="expression", ylabel="TPM")
plt.suptitle("")
plt.savefig(snakemake.output["boxplot_norm"], dpi=300, bbox_inches="tight")
