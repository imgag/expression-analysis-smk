renv::load(snakemake@params[["renv"]])
import::from("magrittr", "%>%")


counts <-
  data.table::fread(snakemake@input[["counts"]])

metadata <-
  readr::read_tsv(snakemake@input[["metadata"]], show_col_types = FALSE)


fpkm_to_tpm <-
  function(fpkm_tbl) {
    prop.table(fpkm_tbl, 2) * 1e6
  }


d <-
  edgeR::DGEList(
    counts = counts[, -c(1, 2)],
    group = metadata$group,
    genes = counts[, 1:2]
  )

# cpm filtered and TMM normalized expression data set
keep <-
  edgeR::filterByExpr(d,
    group = metadata$group,
    min.count = snakemake@params[["min_count"]],
    min.total.count = snakemake@params[["min_total_count"]]
  )

# excluded_genes <-
#   unlist(counts[, 1]) %in% snakemake@config[["exclude_genes"]]

excl_genes <-
  tryCatch(as.character(unlist(read.table(snakemake@input[["excluded_genes"]]))), error=function(e) c())

excluded_genes <-
  unlist(counts[, 1]) %in% excl_genes

d_filtered <-
  d[keep & !excluded_genes, , keep.lib.sizes = FALSE] %>%
  edgeR::calcNormFactors()

# library sizes
library_sizes <-
  tibble::tibble(
    sample = rownames(d$samples),
    full = d$samples$lib.size,
    filtered = d_filtered$samples$lib.size,
    dropped = full - filtered
  )


# output tables
tables_full <-
  list(
    `raw` = cbind(d$genes, d$counts),
    `cpm` = cbind(d$genes, edgeR::cpm(d)),
    `fpkm` = cbind(d$genes, edgeR::rpkm(d)),
    `tpm` = cbind(d$genes, fpkm_to_tpm(edgeR::rpkm(d)))
  )
tables_filtered <-
  list(
    `raw` = cbind(d_filtered$genes, d_filtered$counts),
    `cpm` = cbind(d_filtered$genes, edgeR::cpm(d_filtered)),
    `fpkm` = cbind(d_filtered$genes, edgeR::rpkm(d_filtered)),
    `tpm` = cbind(d_filtered$genes, fpkm_to_tpm(edgeR::rpkm(d_filtered)))
  )

purrr::iwalk(
  tables_full,
  ~ readr::write_tsv(.x, snakemake@output[[sprintf("full_%s", .y)]])
)
purrr::iwalk(
  tables_filtered,
  ~ readr::write_tsv(.x, snakemake@output[[sprintf("filtered_%s", .y)]])
)
readr::write_tsv(library_sizes, snakemake@output[["library_sizes"]])
save(list = c("d_filtered", "d"), file = snakemake@output[["edger_rdata"]])