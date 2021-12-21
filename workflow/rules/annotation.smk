ruleorder: manual_annotation > prepare_annotation


rule prepare_annotation:
    output:
        "intermediate/annotation.tsv",
    params:
        build=config["build"],
        renv=srcdir("../renv"),
    conda:
        "../conda/renv.yml"
    script:
        "../scripts/export_annotables.R"


rule manual_annotation:
    input:
        config["annotation"],
    output:
        "intermediate/annotation.tsv",
    shell:
        """
        cp {input} {output}
        """
