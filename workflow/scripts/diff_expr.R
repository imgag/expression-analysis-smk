renv::load(snakemake@params[["renv"]])
import::from("magrittr", "%>%")

load(snakemake@input[["edger_rdata"]])

metadata <-
  readr::read_tsv(snakemake@input[["metadata"]], show_col_types = FALSE)

contrast_table <-
  readr::read_tsv(snakemake@input[["contrast_table"]], show_col_types = FALSE)

for (fac in names(snakemake@config[["factors"]])) {
  metadata[fac] <-
    factor(metadata[[fac]], levels = snakemake@config[["factors"]][[fac]])
}

design <-
  model.matrix(as.formula(snakemake@params[["model"]]), data = metadata)

colnames(design) <-
  make.names(colnames(design))

print(design)

# estimate dispersions
d_filtered <-
  edgeR::estimateDisp(d_filtered,
    design = design,
    robust = TRUE
  )

# fit generalized linear model
fit <-
  edgeR::glmQLFit(d_filtered, design)


contrasts <-
  limma::makeContrasts(
    contrasts = contrast_table$coefficients,
    levels = design
  )

nb_tests <-
  colnames(contrasts) %>%
  rlang::set_names(contrast_table[match(., contrast_table$coefficients), ]$name) %>%
  purrr::map(~ edgeR::glmQLFTest(fit, contrast = contrasts[, .x]) %>%
    edgeR::topTags(.,
      n = nrow(.),
      adjust.method = "fdr",
      sort.by = "PValue"
    ) %>%
    as.data.frame() %>%
    tibble::as_tibble())

nb_tests %>%
  purrr::iwalk(~ readr::write_tsv(
    .x,
    sprintf("results/diff_expr/contrast.%s.tsv", .y)
  ))

nb_tests_summary <-
  nb_tests %>%
  purrr::imap_dfr(~ dplyr::filter(
    .x,
    abs(logFC) >= snakemake@params[["logfc"]],
    PValue < snakemake@params[["pvalue"]],
    FDR < snakemake@params[["fdr"]]
  ) %>%
    dplyr::mutate(
      direction = factor(
        logFC >= 0,
        levels = c(TRUE, FALSE),
        labels = c("up", "down")
      )
    ) %>%
    dplyr::count(direction, .drop = FALSE),
  .id = "name"
  ) %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  dplyr::mutate(total = up + down)

readr::write_tsv(nb_tests_summary, snakemake@output[["test_summary"]])

save(
  list = c("d_filtered", "fit", "nb_tests"),
  file = snakemake@output[["edger_rdata"]]
)