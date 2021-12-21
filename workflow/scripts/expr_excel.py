from excel_helper import *

wb = make_workbook()

ws = add_worksheet(wb, "cpm", "Gene expression, filtered, CPM")
add_table(
    ws,
    pd.read_csv(snakemake.input["filtered_cpm"], sep="\t"),
    "filtered_cpm",
    start_row=5,
)

ws = add_worksheet(wb, "tpm", "Gene expression, filtered, TPM")
add_table(
    ws,
    pd.read_csv(snakemake.input["filtered_tpm"], sep="\t"),
    "filtered_tpm",
    start_row=5,
)

ws = add_worksheet(wb, "raw", "Gene expression, raw, counts")
add_table(
    ws,
    pd.read_csv(snakemake.input["full_raw"], sep="\t"),
    "full_raw",
    start_row=5,
)

wb.save(snakemake.output[0])
