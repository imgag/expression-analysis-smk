renv::load(snakemake@params[["renv"]])


annot <-
    get(data(
        list = snakemake@params[["build"]],
        package = "annotables",
        envir = environment()
    ))

write.table(
    annot,
    file = snakemake@output[[1]], quote = FALSE, sep = "\t", row.names = FALSE
)