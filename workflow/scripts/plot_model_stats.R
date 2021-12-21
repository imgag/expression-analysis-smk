renv::load(snakemake@params[["renv"]])


load(snakemake@input[["edger_rdata"]])


bcv <-
  sqrt(d_filtered$common.dispersion)
cat(bcv, file = snakemake@output[["bcv_value"]])

png(snakemake@output[["bcv_plot"]])
edgeR::plotBCV(d_filtered)
dev.off()


png(snakemake@output[["disp_plot"]])
edgeR::plotQLDisp(fit)
dev.off()