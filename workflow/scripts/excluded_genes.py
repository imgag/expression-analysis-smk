import pandas as pd


def parse_filter_spec(d, inverse=False):
    if not d:
        return set()

    selected_genes = set()
    filter = []
    for el in d:
        if isinstance(el, str):
            # assume symbol
            filter.append({"symbol": el})
        elif isinstance(el, dict):
            filter.append(el)

    for f in filter:
        assert all([k in annot.columns for k in f.keys()])
        filtered = annot
        for k, v in f.items():
            if not inverse:
                filtered = filtered[filtered[k] == str(v)]
            else:
                filtered = filtered[filtered[k] != str(v)]
        selected_genes.update(list(filtered["gene_id"].unique()))

    return selected_genes


annot = pd.read_csv(snakemake.input["annotation"], sep="\t").rename(
    columns={snakemake.config["annotation_gene_id"]: "gene_id"}
)
excluded_genes = parse_filter_spec(
    snakemake.config["exclude_genes"]
) | parse_filter_spec(snakemake.config["exclude_genes_not_matching"], inverse=True)

with open(snakemake.output["excluded_genes"], "w") as out:
    for g in excluded_genes:
        out.write(f"{g}\n")
