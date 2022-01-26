rule diff_expr:
    input:
        edger_rdata=rules.expr_tables.output.edger_rdata,
        metadata=config["samples"],
        contrast_table=config["contrasts"],
    output:
        expand("results/diff_expr/contrast.{n}.tsv", n=contrasts),
        test_summary="results/diff_expr/summary.tsv",
        edger_rdata="intermediate/edger_diff_expr.rda",
    params:
        model=config["model"],
        pvalue=config["de_filter"]["pvalue"],
        fdr=config["de_filter"]["fdr"],
        logfc=config["de_filter"]["logfc"],
        renv=srcdir("../renv"),
    conda:
        "../conda/renv.yml"
    script:
        "../scripts/diff_expr.R"


rule plot_model_stats:
    input:
        edger_rdata="intermediate/edger_diff_expr.rda",
    output:
        bcv_value="results/diff_expr/bcv.txt",
        bcv_plot="results/diff_expr/plot_bcv.png",
        disp_plot="results/diff_expr/plot_qldisp.png",
    params:
        renv=srcdir("../renv"),
    conda:
        "../conda/renv.yml"
    script:
        "../scripts/plot_model_stats.R"


rule diff_expr_annotated:
    input:
        tbl="results/diff_expr/contrast.{contrast}.tsv",
        annotation="intermediate/annotation.tsv",
    output:
        tbl="results/diff_expr_annot/contrast.{contrast}.tsv",
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


rule diff_expr_volcano:
    input:
        metadata=config["samples"],
        contrasts=config["contrasts"],
        contrast="results/diff_expr_annot/contrast.{contrast}.tsv",
    output:
        plot="results/diff_expr_annot/contrast.{contrast}.volcano.png",
    conda:
        "../conda/plot_sns.yml"
    script:
        "../scripts/volcano.py"


rule diff_expr_heatmap:
    input:
        metadata=config["samples"],
        contrasts=config["contrasts"],
        contrast="results/diff_expr_annot/contrast.{contrast}.tsv",
        filtered_tpm="results/expr/filtered_tpm.tsv",
        group_colors="intermediate/group_colors.tsv",
    output:
        heatmap="results/diff_expr_annot/contrast.{contrast}.heatmap.png",
    conda:
        "../conda/plot_sns.yml"
    script:
        "../scripts/heatmap.py"


rule diff_expr_excel:
    input:
        contrasts=expand("results/diff_expr_annot/contrast.{n}.tsv", n=contrasts),
        heatmaps=expand("results/diff_expr_annot/contrast.{n}.heatmap.png", n=contrasts),
        volcanos=expand("results/diff_expr_annot/contrast.{n}.volcano.png", n=contrasts),
        metadata=config["samples"],
        contrast_table=config["contrasts"],
        test_summary="results/diff_expr/summary.tsv",
        corr_tbl="results/expr/corr_tbl.tsv",
        clustermap="results/expr/clustermap.png",
        pca="results/expr/pca.png",
        library_sizes=rules.expr_tables.output.library_sizes,
        filtered_raw=rules.expr_tables.output.filtered_raw,
    output:
        "results/diff_expr.xlsx",
    conda:
        "../conda/excel.yml"
    script:
        "../scripts/diff_expr_excel.py"
