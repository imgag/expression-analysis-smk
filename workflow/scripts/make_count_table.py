import pandas as pd


def read_featurecounts(path, name):
    df = pd.read_csv(path, sep="\t", comment="#", usecols=[0, 5, 6])
    df.columns = ["gene_id", "gene_length", name]
    df.set_index(["gene_id", "gene_length"], inplace=True)
    return df


sampletable = pd.read_csv(snakemake.input[0], sep="\t")
counts = pd.concat(
    (
        read_featurecounts(p, n)
        for p, n in zip(sampletable["file"], sampletable["sample"])
    ),
    axis=1,
)

counts.reset_index().to_csv(snakemake.output["counts"], sep="\t", index=False)
