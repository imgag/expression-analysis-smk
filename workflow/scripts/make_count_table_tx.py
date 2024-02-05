import pandas as pd


def read_salmon(path, name):
    df = pd.read_csv(path, sep="\t", usecols=[0, 1, 4])
    df.columns = ["transcript_id", "transcript_length", name]
    df.set_index(["transcript_id", "transcript_length"], inplace=True)
    return df


sampletable = pd.read_csv(snakemake.input[0], sep="\t")
counts = pd.concat(
    (
        read_salmon(p, n)
        for p, n in zip(sampletable["file"], sampletable["sample"])
    ),
    axis=1,
)

counts.reset_index().to_csv(snakemake.output["counts"], sep="\t", index=False)
