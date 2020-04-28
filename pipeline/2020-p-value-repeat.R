## All final.tsv are same columns
## Just pick up one to modify it
## The 9-th column is tfbs
#library(data.table)
#dt = fread("/home/zhangjing/roadmap/lung/bed/cal_model_lung_final_result.tsv")
#dt_loc = dt[, .(V1, V2, V3)]
#fwrite(dt_loc, file = "logistic_input_loc.bed", sep = "\t", col.names = FALSE)
## Get the distance to TFBS midpoint
#system("bedtools closest -a logistic_input_loc.bed -b sorted.tfbs_midpoint_3cols.bed -d > dist2tfbs_for_logistic_input.bed")
#
#dist_dt = fread("dist2tfbs_for_logistic_input.bed", header = FALSE)
#dist_dt = dist_dt[, .(V1, V3, V7)]
#str(dist_dt)
#summary(dist_dt$V7)
#
# Assign the max distance to 100
#dist_dt$V7[dist_dt$V7 > 100] = 100
# Assign score to 0-100
#dist_dt$score = 100 - dist_dt$V7
#dist_dt = unique(dist_dt)
#dist_dt$V7 = NULL
#
#dt_update = merge(dt, dist_dt, 
#                  by = c("V1", "V3"), all.x = TRUE)

## Correlation between distance to TFBS midpoint and chip-seq signal intensity
# cor.test(dt_update$score, dt_update$V9)
#dt_update$V9 = dt_update$score
#dt_update$score = NULL

#fwrite(dt_update, file = "logistic_input_2020.tsv", col.names = TRUE, sep = "\t")


modify_tfbs = function(proj_file) {
  library(data.table)
  message("Modifying tfbs column in ", proj_file, "...")
  dt = fread(proj_file)
  dist_dt = fread("dist2tfbs_for_logistic_input.bed", header = FALSE)
  dist_dt = dist_dt[, .(V1, V3, V7)]
  
  # Assign the max distance to 100
  dist_dt$V7[dist_dt$V7 > 100] = 100
  # Assign score to 0-100
  dist_dt$score = 100 - dist_dt$V7
  dist_dt = unique(dist_dt)
  dist_dt$V7 = NULL
  
  dt_update = merge(dt, dist_dt, 
                    by = c("V1", "V3"), all.x = TRUE)
  dt_update$V9 = dt_update$score
  dt_update$score = NULL
  ## Correct the order
  dt_update = dt_update[, c(1, 3, 2, 4:ncol(dt_update)), with = F]
  message("Done. Outputting result file.")
  result_file = paste0("update_", basename(proj_file))
  fwrite(dt_update, 
         file = result_file, 
         col.names = TRUE, sep = "\t")
  gc()
  return(result_file)
}


# Calling logistical models -----------------------------------------------

library(stringr)
library(pROC)
library(caret)
library(data.table)
source("logit_function20181112.R")

type_list = list(
  lung = c("LUSC-CN", "LUSC-KR", "mock"),
  esophagus = c("ESAD-UK","ESCA-CN","mock"),
  liver = c("LIAD-FR", "LICA-CN", "LICA-FR", "LINC-JP", "LIRI-JP","mock"),
  breast = c("BRCA-EU","BRCA-FR","BRCA-US","mock"),
  pancreas = c("PACA-AU","PACA-CA","PAEN-AU","PAEN-IT","mock"),
  kidney = c("RECA-EU","mock"),
  blood = c("ALL-US","CLLE-ES","MALY-DE","NKTL-SG","mock"),
  ovary = c("OV-AU","mock"),
  melanoma = c("MELA-AU","SKCA-BR","SKCM-US","mock")
)

file_list = list(
  v1 = "/home/zhangjing/roadmap/lung/bed/cal_model_lung_final_result.tsv",
  v2 = "/home/zhangjing/roadmap/Esophagus/bed/cal_model_esophagus_final_result.tsv",
  v3 = "/home/zhangjing/roadmap/liver/bed/cal_model_liver_final_result.tsv",
  v4 = "/home/zhangjing/roadmap/breast/bed/cal_model_breast_final_result.tsv",
  v5 = "/home/zhangjing/roadmap/Pancreas/bed/cal_model_pancreas_final_result.tsv",
  v6 = "/home/zhangjing/roadmap/Kidney/bed/cal_model_kidney_final_result.tsv",
  v7 = "/home/zhangjing/roadmap/blood/bed/cal_model_blood_final_result.tsv",
  v8 = "/home/zhangjing/roadmap/Ovary/bed/cal_model_ovary_final_result.tsv",
  v9 = "/home/zhangjing/roadmap/E061_Melanoma/bed/cal_model_mela_final_result.tsv"
)
names(file_list) = names(type_list)

call_model = function(proj_index) {
  timer <- Sys.time()
  message("Calling model for cancer type: ", proj_index)
  result_file = paste0("logit_", proj_index, ".RData")

  #if (file.exists(result_file)) {
 #   message("The result has been called!")
 #   return(NULL)
 # }

  modified_file = modify_tfbs(file_list[[proj_index]])
  mut <- logit_form(type_list[[proj_index]], modified_file)
  if (proj_index %in% c("breast", "pancreas", "kidney", "ovary", "melanoma")) {
    ## Add dnase data
    dnase_file = switch(
      proj_index,
      breast = "/home/zhangjing/roadmap/breast/bed/E028-DNase.hotspot.broad.bed",
      pancreas = "/home/zhangjing/roadmap/Pancreas/bed/E098-DNase.hotspot.broad.bed",
      kidney = "/home/zhangjing/roadmap/Kidney/bed/E086-DNase.hotspot.broad.bed",
      ovary = "/home/zhangjing/roadmap/Ovary/bed/E097-DNase.hotspot.broad.bed",
      melanoma = "/home/zhangjing/roadmap/E061_Melanoma/bed/E059-DNase.hotspot.broad.bed"
    )
    
    message("Adding dnase data...")
    mut <- add_dnase(mut, dnase_file)
  }
  loc_info = mut[, c("donor", "chr", "end")]
  message("print loc...")
  print(head(loc_info))
  print(dim(loc_info))
  save(loc_info, file = paste0("logit_loc_", proj_index, ".RData"))
  mut <- format_trans(mut)
  message("print model input data...")
  print(head(mut))
  print(dim(mut))
  save(mut, file = paste0("logit_input_", proj_index, ".RData"))
  file.remove(modified_file)
  #cross_val(mut, result_file)
  if (file.exists(result_file)) {
    message("The result has been called!")
    return(NULL)
  }
  cross_val(mut, result_file)
  message("Done")
  gc()
  #file.remove(modified_file)
  message(Sys.time() - timer)
}

cancer_type = c("blood", "breast", "esophagus",
                "kidney", "liver", "lung",
                "ovary", "pancreas", "melanoma")

for (i in cancer_type) {
  call_model(i)
}
