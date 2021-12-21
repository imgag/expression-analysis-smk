from excel_helper import *

wb = make_workbook()


# Sample Overview
ws = add_worksheet(wb, "samples", "Sample Overview")
add_table(
    ws, pd.read_csv(snakemake.input["metadata"], sep="\t").drop(columns=["file"]), "metadata", start_row=5
)


# Sample Similarity
ws = add_worksheet(wb, "similarity", "Sample Similarity")
add_table(
    ws, pd.read_csv(snakemake.input["corr_tbl"], sep="\t"), "corr_tbl", start_row=5
)
add_image(ws, f"A{ws._current_row + 2}", snakemake.input["clustermap"], width=800)
add_image(ws, f"G{ws._current_row + 2}", snakemake.input["pca"], width=800)


# Expression Filtering
ws = add_worksheet(wb, "filter", "Expression Filter")
ws._current_row += 1
ws.append(
    [
        "Minimum count in worthwhile number of samples:",
        snakemake.config["cpm_filter"]["min_count"],
    ]
)
ws.append(["Minimum total count:", snakemake.config["cpm_filter"]["min_total_count"]])
ws.append(
    [
        "Number of retained genes:",
        len(pd.read_csv(snakemake.input["filtered_raw"], sep="\t").index),
    ]
)
add_table(
    ws,
    pd.read_csv(snakemake.input["library_sizes"], sep="\t"),
    "library_sizes",
    start_row=7,
    col_type="@iii",
)


# Contrast Overview
ws = add_worksheet(wb, "contrasts", "Contrast Overview")
add_table(
    ws,
    pd.read_csv(snakemake.input["contrast_table"], sep="\t"),
    "contrast_table",
    start_row=5,
)


# Differential Expression Summary
ws = add_worksheet(wb, "summary", "Differential Expression Summary")
ws._current_row += 1
ws.append(["Model:", snakemake.config["model"]])
add_table(
    ws,
    pd.read_csv(snakemake.input["test_summary"], sep="\t"),
    "test_summary",
    start_row=5,
)


# Sheet for each contrast
contrasts = pd.read_csv(snakemake.input["contrast_table"], sep="\t", index_col=0)

icons_fdr = IconSetRule(
    "3Symbols",
    "num",
    [0, snakemake.config["de_filter"]["fdr"], snakemake.config["de_filter"]["fdr"]],
    reverse=True,
)

for n, row in contrasts.iterrows():
    ws = add_worksheet(wb, n, row["title"])
    ws["A3"] = row["description"]
    df = pd.read_csv(f"results/diff_expr_annot/contrast.{n}.tsv", sep="\t")
    columns = (
        ["gene_id"]
        + snakemake.config["annotation_columns"][1:]
        + ["logFC", "PValue", "FDR"]
    )
    add_table(
        ws,
        df.loc[:, columns],  # todo
        f"contrast_{n}",
        start_row=5,
        col_type="@" * (len(columns) - 3) + "nee",
    )
    ws.conditional_formatting.add(
        make_range("E6", (len(df.index), 1)), logfc_colorscale
    )
    ws.conditional_formatting.add(make_range("G6", (len(df.index), 1)), icons_fdr)

    add_image(ws, "I5", f"results/diff_expr_annot/contrast.{n}.volcano.png", height=300)
    add_image(ws, "I21", f"results/diff_expr_annot/contrast.{n}.heatmap.png", width=500)


for ws in wb:
    set_column_width(ws, 20)

wb.save(snakemake.output[0])
