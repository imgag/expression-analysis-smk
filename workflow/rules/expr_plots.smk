rule group_colormap:
    input:
        metadata=config["samples"],
    output:
        colormap=temporary("intermediate/group_colors.tsv"),
    conda:
        "../conda/plot.yml"
    script:
        "../scripts/colormap.py"


rule expr_dist:
    input:
        filtered_cpm="results/expr/filtered_cpm.tsv",
        full_cpm="results/expr/full_cpm.tsv",
    output:
        boxplot="results/expr/cpm_boxplot.png",
        density="results/expr/cpm_density.png",
    conda:
        "../conda/plot.yml"
    script:
        "../scripts/expr_dist.py"


rule expr_corr:
    input:
        filtered_raw="results/expr/filtered_raw.tsv",
    output:
        table="results/expr/corr_tbl.tsv",
        dendrogram="results/expr/corr_dendrogram.png",
    conda:
        "../conda/plot.yml"
    script:
        "../scripts/expr_corr.py"


rule expr_clustermap:
    input:
        metadata=config["samples"],
        corr_tbl="results/expr/corr_tbl.tsv",
        group_colors="intermediate/group_colors.tsv",
    output:
        clustermap="results/expr/clustermap.png",
    conda:
        "../conda/plot_sns.yml"
    script:
        "../scripts/expr_clustermap.py"


rule expr_pca:
    input:
        metadata=config["samples"],
        filtered_cpm="results/expr/filtered_cpm.tsv",
        group_colors="intermediate/group_colors.tsv",
    output:
        pca_plt="results/expr/pca.png",
        pca_tbl="results/expr/pca.tsv",
    conda:
        "../conda/plot_sns_pca.yml"
    script:
        "../scripts/pca.py"


rule expr_gene_boxplots:
    input:
        metadata=config["samples"],
        filtered_tpm_annot="results/expr_annot/filtered_tpm.tsv",
    output:
        boxplot_log="results/genes/{gene}.boxplot.log2.png",
        boxplot_norm="results/genes/{gene}.boxplot.png",
    conda:
        "../conda/plot_sns.yml"
    script:
        "../scripts/boxplot.py"
