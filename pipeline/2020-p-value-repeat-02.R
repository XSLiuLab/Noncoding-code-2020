library(data.table)

cancer_type = c("blood", "breast", "esophagus",
                "kidney", "liver", "lung",
                "ovary", "pancreas", "melanoma")

for (i in cancer_type) {
  message("Processing ", i, "...")
  file_names <- list.files(pattern = paste0(i, ".RData"))
  for (fl in file_names) {
    load(fl)
  }
  loc_info$pred <- as.numeric(predict.glm(prei, type = "response", newdata = mut))
  rm(prei, mut); gc()
  fwrite(loc_info, file = paste0("pred_", i, ".tsv.gz"), sep = "\t")
  rm(loc_info); gc()
}
