
rule make_count_table:
    input:
        config["samples"],
    output:
        counts="intermediate/counts.tsv",
    script:
        "../scripts/make_count_table.py"

rule make_count_table_tx:
    input:
        config["samples"],
    output:
        counts="intermediate/counts.tsv",
    script:
        "../scripts/make_count_table_tx.py"


rule expr_tables:
    input:
        counts="intermediate/counts.tsv",
        metadata=config["samples"],
        excluded_genes="intermediate/excluded_genes.txt",
    output:
        library_sizes="results/expr/library_sizes.tsv",
        full_raw="results/expr/full_raw.tsv",
        full_cpm="results/expr/full_cpm.tsv",
        full_fpkm="results/expr/full_fpkm.tsv",
        full_tpm="results/expr/full_tpm.tsv",
        filtered_raw="results/expr/filtered_raw.tsv",
        filtered_cpm="results/expr/filtered_cpm.tsv",
        filtered_fpkm="results/expr/filtered_fpkm.tsv",
        filtered_tpm="results/expr/filtered_tpm.tsv",
        edger_rdata="intermediate/edger_expr.rda",
    conda:
        "../conda/renv.yml"
    params:
        gene_id=config["annotation_gene_id"],
        min_count=config["cpm_filter"]["min_count"],
        min_total_count=config["cpm_filter"]["min_total_count"],
        renv=srcdir("../renv"),
    script:
        "../scripts/expr_tables.R"


rule expr_tables_annotated:
    input:
        tbl="results/expr/{type}.tsv",
        annotation="intermediate/annotation.tsv",
    output:
        tbl="results/expr_annot/{type}.tsv",
    params:
        gene_id=config["annotation_gene_id"],
        annotation_columns=config["annotation_columns"],
    run:
        annot = pd.read_csv(input["annotation"], sep="\t")
        tbl = pd.read_csv(input["tbl"], sep="\t")
        tbl_annot = annotate(
            tbl, annot, "gene_id", params["gene_id"], params["annotation_columns"]
        )
        tbl_annot.to_csv(output["tbl"], sep="\t", index=False)


rule expr_excel:
    input:
        filtered_cpm="results/expr_annot/filtered_cpm.tsv",
        filtered_tpm="results/expr_annot/filtered_tpm.tsv",
        full_raw="results/expr_annot/full_raw.tsv",
    output:
        "results/expression.xlsx",
    conda:
        "../conda/excel.yml"
    script:
        "../scripts/expr_excel.py"


rule make_excluded_genes:
    input:
        annotation="intermediate/annotation.tsv",
    output:
        excluded_genes="intermediate/excluded_genes.txt",
    script:
        "../scripts/excluded_genes.py"
