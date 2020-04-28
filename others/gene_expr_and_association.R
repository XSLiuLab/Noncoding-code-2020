library(UCSCXenaTools)
library(dplyr)
# ETS transcripts and CDC20
ETS_list = c("ERF", "ETV3", "ELF3", "ELF5", "ETV4",
             "ETV5", "ETV1", "ELK1", "ELK4", "ELK3",
             "ETV6", "ETV7")
gene_list = c(ETS_list, "CDC20")

cohort = XenaData %>%
  dplyr::filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan("TCGA Melanoma")   # select LUAD cohort
cohort

cli_query = cohort %>%
  dplyr::filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload()


table(cli$sample_type)
cli = XenaPrepare(cli_query)
tumor_id = cli$sampleID[cli$sample_type %in% c("Additional Metastatic", "Metastatic", "Primary Tumor")]
meta_id = cli$sampleID[cli$sample_type %in% c("Additional Metastatic", "Metastatic")]
pri_id = cli$sampleID[cli$sample_type %in% c("Primary Tumor")]

ge = cohort %>%
  dplyr::filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq")

expr_mat = fetch_dense_values(host = ge$XenaHosts,
                          dataset = ge$XenaDatasets,
                          identifiers = gene_list,
                          use_probeMap = TRUE)
# From GDC hub
expr_mat2 = fetch_dense_values(host = "https://gdc.xenahubs.net",
                              dataset = "TCGA-SKCM.htseq_fpkm.tsv",
                              identifiers = gene_list,
                              use_probeMap = TRUE)
colnames(expr_mat2) = substr(colnames(expr_mat2), 1, 15)

rownames(expr_mat)
rownames(expr_mat2)

tumor_expr = expr_mat[, colnames(expr_mat) %in% tumor_id]
meta_expr = expr_mat[, colnames(expr_mat) %in% meta_id]
pri_expr = expr_mat[, colnames(expr_mat) %in% pri_id]

tumor_expr2 = expr_mat2[, colnames(expr_mat2) %in% tumor_id]
meta_expr2 = expr_mat2[, colnames(expr_mat2) %in% meta_id]
pri_expr2 = expr_mat2[, colnames(expr_mat2) %in% pri_id]


plot_boxplot_and_heatmap = function(tumor_expr, fn, title = "Expression and Correlation Heatmap") {
  pdf(fn)
  
  boxplot(t(tumor_expr), las=2, main = title)
  #mtext(title, outer=TRUE,  cex=1)
  cor_mat = cor(t(tumor_expr))
  breaksList <- seq(-1, 1, by = 0.01)
  pheatmap::pheatmap(cor_mat,
                     display_numbers = TRUE,
                     color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                     breaks = breaksList
  )
  
 
  on.exit(dev.off())
}


plot_boxplot_and_heatmap(tumor_expr, 
                         "MELA-ALL-tumors-data-version-2017-10.pdf",
                         title = "Expression and Correlation Heatmap for All Tumors (472 samples)")
plot_boxplot_and_heatmap(meta_expr, 
                         "MELA-Metastatic-tumors-data-version-2017-10.pdf",
                         title = "Expression and Correlation Heatmap for Metastatic Tumors (369 samples)")
plot_boxplot_and_heatmap(pri_expr, 
                         "MELA-Primary-tumors-data-version-2017-10.pdf",
                         title = "Expression and Correlation Heatmap for Primary Tumors (103 samples)")

plot_boxplot_and_heatmap(tumor_expr2, 
                         "MELA-ALL-tumors-data-version-2019-07.pdf",
                         title = "Expression and Correlation Heatmap for All Tumors (471 samples)")
plot_boxplot_and_heatmap(meta_expr2, 
                         "MELA-Metastatic-tumors-data-version-2019-07.pdf",
                         title = "Expression and Correlation Heatmap for Metastatic Tumors (368 samples)")
plot_boxplot_and_heatmap(pri_expr2, 
                         "MELA-Primary-tumors-data-version-2019-07.pdf",
                         title = "Expression and Correlation Heatmap for Primary Tumors (103 samples)")



