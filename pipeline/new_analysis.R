# New Analysis BEGIN ------------------------------------------------------
# 计算1Mbp区间内突变频率与遗传、表观注释值的相关性
# Coding, Noconding, Promoter ...
# Coding为外显子区域

setwd("~/New_Analysis/")
source("functions.R")
library(data.table)
Sys.setenv(PATH=paste0("/home/zhangjing/bedtools2/bin:", Sys.getenv("PATH")))

# 预处理数据 -------------------------------------------------------------------
options(scipen = 30)
original_mut = "../icgc_data/预处理所有突变.tsv"    # 经过预处理过的突变记录文件

if (F) {
  mut <- fread(original_mut, data.table = F)
  pos_mut <- mut[mut$mut == "single base substitution",]
  ###以bed格式输出所有单碱基突变
  mut_bed <- pos_mut[,3:5]
  mut_bed$start <- mut_bed$start - 1L
  mut_bed$start = as.integer(mut_bed$start)
  mut_bed$end = as.integer(mut_bed$end)
  write.table(mut_bed, "single_mut.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  rm(mut, pos_mut, mut_bed); gc()
}

# 编码区 ---------------------------------------------------------------------

get_cds_region("../predict_prob/Homo_sapiens.GRCh37.75.gtf")
gc()
# merge，以防重叠区间
system("sort -k 1,1 -k 2,2n cds_region.bed > sorted_cds_region.bed")
#Sys.getenv("PATH")
system("bedtools merge -i sorted_cds_region.bed > merged_sorted_cds_region.bed")
file.remove("sorted_cds_region.bed")

# 非编码区 --------------------------------------------------------------------
system("sort -k 1,1 -k 2,2n human.hg19.genome > sorted_human.hg19.genome")  #bedtools包中附带文件
system("bedtools complement -i merged_sorted_cds_region.bed -g sorted_human.hg19.genome > region_noncoding.bed")

# Promoter 区 --------------------------------------------------------------
# 选取-2500 到 +500 作为启动子区域
up_site = 2500
down_site = 500
##过滤数据
gtf_ensmble <- fread("../predict_prob/Homo_sapiens.GRCh37.75.gtf", 
                     sep = "\t", colClasses=list(character= 1))
names(gtf_ensmble)[1] <- "V1"
ensmble_pro_cod <- gtf_ensmble[gtf_ensmble$V2 == "protein_coding",] #选择仅仅是蛋白质编码的转录本
transcript_pro_cod_ensmble <- ensmble_pro_cod[ensmble_pro_cod$V3 == "transcript",]          ##########从转录本中选取总的转录本区间
options(scipen=30)
chr_index <- c(c(1:22),"X","Y")
chr_trans_pro <- transcript_pro_cod_ensmble[transcript_pro_cod_ensmble$V1 %in% chr_index, ] #仅仅挑选chr1:22 + X + Y 染色体
##从SST获取up_site到down_site的区间
pos_trans_chr_gtf <- chr_trans_pro[chr_trans_pro$V7 == "+",]
neg_trans_chr_gtf <- chr_trans_pro[chr_trans_pro$V7 == "-",]
pos_chr <- paste("chr",pos_trans_chr_gtf$V1, sep = "")
pos_start <- pos_trans_chr_gtf$V4 - up_site
pos_end <- pos_trans_chr_gtf$V4 + down_site
pos_bed <- data.frame(pos_chr, pos_start, pos_end)
names(pos_bed) <- c("chr","start","end")

neg_chr <- paste("chr",neg_trans_chr_gtf$V1, sep = "")
neg_start <- neg_trans_chr_gtf$V5 + up_site
neg_end <- neg_trans_chr_gtf$V5 - down_site
neg_bed <- data.frame(neg_chr, neg_end, neg_start)
names(neg_bed) <- c("chr","start","end")

promoter_region <- rbind(pos_bed, neg_bed)
promoter_region$index <- 1
write.table(promoter_region,"promoter_region.bed", sep = "\t", row.names = F, col.names = F,quote = F)

# 同样地，排序并merge
system("sort -k 1,1 -k 2,2n promoter_region.bed > sorted_promoter_region.bed")
system("bedtools merge -i sorted_promoter_region.bed > merged_sorted_promoter_region.bed")
file.remove("sorted_promoter_region.bed")

# 非Promoter的noncoding区 ----------------------------------------------------
# noncoding区间除了Promoter之外的区间

system("bedtools subtract -a region_noncoding.bed -b merged_sorted_promoter_region.bed > region_non_promoter.bed")


rm(list = ls()); gc()
source("functions.R")
# 获取不同区间上的突变数目 --------------------------------------------------------
region.coding = "merged_sorted_cds_region.bed"
region.noncoding = "region_noncoding.bed"
region.promoter = "merged_sorted_promoter_region.bed"
region.nonpromoter = "region_non_promoter.bed"

original_bed = "single_mut.bed"

get_reg_mution(original_bed, region.coding, "mut_coding.bed")
get_reg_mution(original_bed, region.noncoding, "mut_noncoding.bed")
get_reg_mution(original_bed, region.promoter, "mut_promoter.bed")
get_reg_mution(original_bed, region.nonpromoter, "mut_nonpromoter.bed")


# 分割基因组区间 -----------------------------------------------------------------
get_region_bed(1000000L)
system("sort -k 1,1 -k 2,2n hg19_mb_1M.bed > sorted_hg19_1M.bed")

# 计算突变数 -------------------------------------------------------------------
get_reg_mution_number("mut_coding.bed", "sorted_hg19_1M.bed", "1M_mut_coding.tsv")
get_reg_mution_number("mut_noncoding.bed", "sorted_hg19_1M.bed", "1M_mut_noncoding.tsv")
get_reg_mution_number("mut_promoter.bed", "sorted_hg19_1M.bed", "1M_mut_promoter.tsv")
get_reg_mution_number("mut_nonpromoter.bed", "sorted_hg19_1M.bed", "1M_mut_nonpromoter.tsv")


# 计算突变的各种注释结果 -------------------------------------------------------------

# 这里取自师兄之前的结果
# 下载原数据
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/mutation_counts_mb.tsv",
              destfile = "mutation_counts_mb.tsv")
# CpG
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/CpG_percent_Mb.bed",
              destfile = "CpG_percent_Mb.bed")
# GC
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/hg19_gcMb.bed",
              destfile = "hg19_gcMb.bed")

# Conservation
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/conservation_Mb.bed",
              destfile = "conservation_Mb.bed")

# Reptime
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/reptime_Mb.bed",
              destfile = "reptime_Mb.bed")

# polII
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/polII_mb.bed",
              destfile = "polII_mb.bed")

# mappability
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/mean_mb_mappability.tsv",
              destfile = "mean_mb_mappability.tsv")

# rec_rate
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/rec_rate.tsv",
              destfile = "rec_rate.tsv")

# tfbs
download.file("https://raw.githubusercontent.com/zhangjing1589/noncoding/master/tfbs_mb.bed",
              destfile = "tfbs_mb.bed")

## Added on 2020
#system("sort -k 1,1 -k 2,2n tfbs_midpoint.union.bed > sorted.tfbs_midpoint.bed")

dt = data.table::fread("tfbs_midpoint.union.bed", header = FALSE)
dt$V2 = as.integer(dt$V2)
dt$V3 = as.integer(dt$V3)
data.table::fwrite(unique(dt[, c(1, 2, 3)]), 
                   file = "tfbs_midpoint_3cols.bed",
                   sep = "\t",
                   col.names = FALSE)
system("sort -k 1,1 -k 2,2n tfbs_midpoint_3cols.bed > sorted.tfbs_midpoint_3cols.bed")
system("bedtools closest -a mut_noncoding.bed -b sorted.tfbs_midpoint_3cols.bed -d > dist2tfbs_midpoint.bed")

dt = data.table::fread("dist2tfbs_midpoint.bed", header = FALSE)
summary(dt$V7)
dt = dt[V7 != -1]

dt$V7[dt$V7 > 100] = 100
dt$V4 = NULL
dt$V5 = NULL
dt$V6 = NULL

data.table::fwrite(dt, file = "dist2tfbs.bed", col.names = FALSE, sep = "\t")
rm(list = ls()); gc()

system("bedtools map -a sorted_hg19_1M.bed -b dist2tfbs.bed -c 4 -o count,sum > tfbs_mb_dist.bed")

dist = data.table::fread("tfbs_mb_dist.bed", header = FALSE)

str(dist)

dist$V6 = as.integer(dist$V6)
dist$V6 = dist$V6 / dist$V5

summary(dist$V6)

cor.test(dist$V5, dist$V6)
plot(dist$V6, dist$V5,
     ylab = "Noncoding mutactions in 1Mb region",
     xlab = "Mean distance (bp) to midpoint of TFBS in 1Mb region")

data.table::fwrite(dist, file = "tfbs_mb_dist_tidy.bed", col.names = FALSE, sep = "\t")

# plot mutation density around tfbs
dt = data.table::fread("dist2tfbs_midpoint.bed", header = FALSE)
summary(dt$V7)
dt = dt[V7 != -1]

n_tfbs = length(unique(paste(dt$V4, dt$V6, sep = "-")))

dt$dist = dt$V3 - dt$V6
dist = dt$dist[dt$V7 <= 100]

save(dist, file = "tfbs_dist.RData")
rm(list = ls()); gc()

load(file = "tfbs_dist.RData")
layout(matrix(1:2, nrow=1))
hist(dist, breaks = 50,
     xlab = "Distance (bp) to midpoint of TFBS", main = "Mutation histogram")
plot(density(dist), xlab = "Distance (bp) to mcidpoint of TFBS", main = "Mutation density")
layout(1)

df = table(dist)
df2 = data.table::data.table(
  x = as.integer(names(df)),
  y = as.integer(df) / n_tfbs
)

save(df2, file = "df2.RData")

library(ggplot2)
library(cowplot)
p = ggplot(as.data.frame(df2), aes(x=x, y=y)) + 
  geom_line() + labs(x = "Position relative to TFBS midpoint (bp)",
                     y = "Mutation density") + cowplot::theme_cowplot()

save_plot("Mutation_density_vs_TFBS_dist.pdf",
          plot = p, base_aspect_ratio = 1.6)  

##


# 收集遗传变异相关性 ---------------------------------------------------------------
rm(list = ls());gc()

mut_coding      = "1M_mut_coding.tsv"
mut_noncoding  = "1M_mut_noncoding.tsv"
# mut_noncoding   = "mutation_counts_mb.tsv" #师兄使用的数据
mut_promoter    = "1M_mut_promoter.tsv"
mut_nonpromoter = "1M_mut_nonpromoter.tsv"


annot_CpG          = "CpG_percent_Mb.bed"
annot_GC           = "hg19_gcMb.bed"
annot_conservation = "conservation_Mb.bed"
annot_Reptime      = "reptime_Mb.bed"
annot_polII        = "polII_mb.bed"
annot_mappability  = "mean_mb_mappability.tsv"
annot_rec_rate     = "rec_rate.tsv"
# Modify TFBS to distance from mutation to TFBS midpoint
annot_tfbs         = "tfbs_mb.bed"

annot = c(annot_CpG, annot_GC, annot_conservation, annot_Reptime,
          annot_polII, annot_mappability, annot_rec_rate, annot_tfbs)

region.coding = "merged_sorted_cds_region.bed"
region.noncoding = "region_noncoding.bed"
region.promoter = "merged_sorted_promoter_region.bed"
region.nonpromoter = "region_non_promoter.bed"

add_coverage = function(gn_region, mut_file, result_file, ref_region="sorted_hg19_1M.bed"){
  tmp_file = paste0("tmp_", basename(gn_region))
  system(paste("bedtools coverage -a", ref_region, "-b", gn_region, ">", tmp_file))
  tmp = fread(tmp_file)
  file.remove(tmp_file)
  mut = fread(mut_file)
  res = merge(mut, tmp, by="V4")
  write.table(res[, c(2:4,1,5,12)], result_file, sep = "\t", row.names = F, col.names = F, quote = F)
}

add_coverage(region.coding, mut_coding, "1M_mut_coding_add_cov.tsv")
add_coverage(region.noncoding, mut_noncoding, "1M_mut_noncoding_add_cov.tsv")
add_coverage(region.promoter, mut_promoter, "1M_mut_promoter_add_cov.tsv")
add_coverage(region.nonpromoter, mut_nonpromoter, "1M_mut_nonpromoter_add_cov.tsv")

# calculate coverage
#system("bedtools coverage -a sorted_hg19_1M.bed -b merged_sorted_cds_region.bed   > test_cov2.bed")

# join by x with V4
# V4, V4, V1, V4, V1, V1, V4, V1

obtain_cor = function(annot_list, join_cols, mut_df) {
  require(data.table)
  out = data.table()
  for (i in seq_along(annot_list)) {
    message("Runing annotation #", i, ": ", annot_list[i])
    value = fread(annot_list[i])
    mut = fread(mut_df)
    colnames(mut)[6] = "cov"
    result = merge(value, mut, by.x = colnames(value)[join_cols[i]], by.y = "V4")

    # calculate correlation
    if (i == 5) {
      result = result[(V3.x != 0) & (V5 !=0)]
      tmp = data.table(
        coeff = cor(result$V3.x, result$V5 / result$cov),
        p_val = cor.test(result$V3.x, result$V5 / result$cov)$p.value
      )
    } else if (i == 7) {
      result = result[(avg != 0) & (V5.y != 0)]
      tmp = data.table(
        coeff = cor(result$avg, result$V5.y / result$cov),
        p_val = cor.test(result$avg, result$V5.y / result$cov)$p.value
      )
    } else {
      result = result[(V5.x != 0) & (V5.x != ".") & (V5.y != 0) & (V5.y != ".") ]
      result[, V5.x:=as.numeric(V5.x)]
      result[, V5.y:=as.numeric(V5.y)]
      tmp = data.table(
        coeff = cor(result$V5.x, result$V5.y / result$cov),
        p_val = cor.test(result$V5.x, result$V5.y / result$cov)$p.value
      )
    }

    out = rbind(out, tmp)
  }
  
  out
}

mut_coding      = "1M_mut_coding_add_cov.tsv"
mut_noncoding  = "1M_mut_noncoding_add_cov.tsv"
# mut_noncoding   = "mutation_counts_mb.tsv" #师兄使用的数据
mut_promoter    = "1M_mut_promoter_add_cov.tsv"
mut_nonpromoter = "1M_mut_nonpromoter_add_cov.tsv"

cor_genetic = list()
cor_genetic[["noncoding"]] = obtain_cor(annot, 
                                        join_cols = c(4, 4, 1, 4, 1, 1, 4, 1), 
                                        mut_noncoding)

cor_genetic[["coding"]] = obtain_cor(annot, 
                                        join_cols = c(4, 4, 1, 4, 1, 1, 4, 1), 
                                        mut_coding)
cor_genetic[["promoter"]] = obtain_cor(annot, 
                                     join_cols = c(4, 4, 1, 4, 1, 1, 4, 1), 
                                     mut_promoter)
cor_genetic[["nonpromoter"]] = obtain_cor(annot, 
                                       join_cols = c(4, 4, 1, 4, 1, 1, 4, 1), 
                                       mut_nonpromoter)


cor_genetic[["noncoding"]]
cor_genetic[["coding"]]
cor_genetic[["promoter"]]
cor_genetic[["nonpromoter"]]

cor_genetic_df = rbindlist(cor_genetic, idcol = TRUE)
setnames(cor_genetic_df, ".id", "region")
cor_genetic_df[, features:=rep(c("CpG island", "GC content",
                                "Conservation", "Replication time",
                                "DNA poly II", "Mappability",
                                "Recombination rate", "TFBS"), 4)]
cor_genetic_df$features = factor(cor_genetic_df$features,
                                 c("Mappability", "Replication time",
                                   "TFBS", "GC content",
                                   "CpG island", "DNA poly II",
                                   "Conservation", 
                                   "Recombination rate"))

## Remove "Noncoding-promoter"
load("pipeline/cor_genetic_df.RData")

library(ggplot2)
library(cowplot)
library(RColorBrewer)

ggplot(cor_genetic_df[features != "TFBS"][region != "nonpromoter"],
       aes(x = features, y = coeff, fill=region)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Paired", 
                    labels = c("Coding", "Noncoding", 
                               "Promoter")) +
  labs(x = "Genetic features", y = "Correlation coefficient", fill = "Region") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p_genetic

p_genetic

save_plot("Genetic_corrplot_rm_TFBS.pdf", plot = p_genetic, base_aspect_ratio = 1.6)  

save(cor_genetic_df, file = "cor_genetic_df.RData")

cor_genetic_df2 = cor_genetic_df
setnames(cor_genetic_df2,
         colnames(cor_genetic_df2),
         c("Genomic region", "Correlation coefficient",
           "P value", "Genetic feature"))
cor_genetic_df2[, `Genomic region` := 
                  ifelse(`Genomic region` == "nonpromoter",
                         "Noncoding-promoter", `Genomic region`)]

write.csv(as.data.frame(cor_genetic_df2), file = "遗传相关性结果3-wsx.csv", quote = FALSE, row.names = FALSE)

# library(gt)
# gt_genetic <- gt(data = cor_genetic_df2)
# gt_genetic 

# value <- fread("CpG_percent_Mb.bed")
# mut <- fread("mutation_counts_mb.tsv")
# result <- merge(value, mut, by.x = "V4", by.y = "V4")
# result <- result[result$V5.x != 0,]
# cor(result$V5.x, result$V5.y)
# cor.test(result$V5.x, result$V5.y)
# 
# 
# mut_num <- read.table("mutation_counts_mb.tsv", stringsAsFactors = F)
# mut_num <- mut_num[mut_num$V5 != ".",]
# mut_num$V5 <- as.numeric(mut_num$V5)
# 
# reptime <- read.table("reptime_Mb.bed", stringsAsFactors = F)
# reptime <- reptime[reptime$V5 != ".",]
# reptime$V5 <- as.numeric(reptime$V5)
# result <- merge(reptime, mut_num, by.x = "V4", by.y = "V4")
# 
# result <- result[result$V5.x != 0, ]
# cor(result$V5.x, result$V5.y)
# cor.test(result$V5.x, result$V5.y)



# 表观相关性预处理 ----------------------------------------------------------------
rm(list = ls()); gc()

###按不同组织将肿瘤分类
mut <- fread("~/icgc_data/预处理所有突变.tsv")
donor <- fread("~/icgc_data/donor.tsv")
index <- donor[,1:2]

merge_result <- merge(mut, index, by.x = "donor", by.y = "icgc_donor_id")

dir.create("epi")

#将不同的肿瘤类型输出
#肝癌
liver <- c("LIAD-FR", "LICA-CN", "LICA-FR", "LINC-JP", "LIRI-JP")
liver_cancer <- merge_result[merge_result$project_code %in% liver,]
write.table(liver_cancer, "epi/liver_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(liver_cancer)

# #胃癌
# gastric <- c("STAD-US")
# gastric_cancer <- merge_result[merge_result$project_code %in% gastric,]
# write.table(gastric_cancer, "gastric_cancer_info.bed", sep = "\t", row.names = F, quote = F)
# rm(gastric_cancer)

#肾癌
kidney <- c("RECA-EU")
kidney_cancer <- merge_result[merge_result$project_code %in% kidney,]
write.table(kidney_cancer, "epi/kidney_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(kidney_cancer)

#肺癌
lung <- c("LUSC-CN","LUSC-KR")
lung_cancer <- merge_result[merge_result$project_code %in% lung,]
write.table(lung_cancer, "epi/lung_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(lung_cancer)

#乳腺癌
breast <- c("BRCA-EU","BRCA-FR","BRCA-US")
breast_cancer <- merge_result[merge_result$project_code %in% breast,]
write.table(breast_cancer, "epi/breast_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(breast_cancer)

#胰腺癌
pancreas <- c("PACA-AU","PACA-CA","PAEN-AU","PAEN-IT")
pancreas_cancer <- merge_result[merge_result$project_code %in% pancreas,]
write.table(pancreas_cancer, "epi/pancreas_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(pancreas_cancer)

#卵巢癌
ovary <- c("OV-AU")
ovary_cancer <- merge_result[merge_result$project_code %in% ovary,]
write.table(ovary_cancer, "epi/ovary_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(ovary_cancer)

#黑色素瘤
mela <- c("MELA-AU","SKCA-BR","SKCM-US")
mela_cancer <- merge_result[merge_result$project_code %in% mela,]
write.table(mela_cancer, "epi/mela_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(mela_cancer)

#食道癌
esophagus <- c("ESAD-UK","ESCA-CN")
esophagus_cancer <- merge_result[merge_result$project_code %in% esophagus,]
write.table(esophagus_cancer, "epi/esophagus_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(esophagus_cancer)

#血癌
blood <- c("ALL-US","CLLE-ES","MALY-DE","NKTL-SG")
blood_cancer <- merge_result[merge_result$project_code %in% blood,]
write.table(blood_cancer, "epi/blood_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(blood_cancer)

# #结肠癌
# colon <- c("COAD-US","COCA-CN")
# colon_cancer <- merge_result[merge_result$project_code %in% colon,]
# write.table(colon_cancer, "colon_cancer_info.bed", sep = "\t", row.names = F, quote = F)
# rm(colon_cancer)
# 
# #胶质瘤
# brain <- c("LGG-US")
# brain_cancer <- merge_result[merge_result$project_code %in% brain,]
# write.table(brain_cancer, "brain_cancer_info.bed", sep = "\t", row.names = F, quote = F)
# rm(brain_cancer)

rm(list = ls()); gc()


#删除chrM染色体
rm_chrM <- function(file)
{
  df <- fread(file, data.table = F)
  df <- df[df$V1 != "chrM",]
  write.table(df, file, sep = "\t", row.names = F, quote = F, col.names = F)
}


##将表观遗传bed文件转化为bigwig文件
bed2bigwig <- function(bedfile, result_file)
{
  temp1 <- paste("sorted_", bedfile, sep = "")
  cmd1 <- paste("sort -k 1,1 -k 2,2n ",bedfile," > ",temp1, sep = "")
  system(cmd1)
  rm_chrM(temp1)
  
  cmd2 <- paste("~/bedGraphToBigWig ",temp1, 
                "~/New_Analysis/sorted_human.hg19.genome",
                file.path("~/New_Analysis/bigwig", result_file))
  system(cmd2)
  file.remove(temp1)
}

# #将DNase文件转换为bed 格式
# toBed <- function(file)
# {
#   df <- fread(file, data.table = F)
#   df <- df[,c(1:3,5)]
#   write.table(df, file, sep = "\t", row.names = F, quote = F, col.names = F)
# }

dir.create("~/New_Analysis/bigwig")

setwd("~/roadmap/lung/bed/")
#肺癌 E096
##将注释文件换成二进制
bed2bigwig("E096-H3K27ac.broadPeak.bed","E096-H3K27ac.broadPeak.bigwig")
bed2bigwig("E096-H3K27me3.broadPeak.bed","E096-H3K27me3.broadPeak.bigwig")
bed2bigwig("E096-H3K36me3.broadPeak.bed","E096-H3K36me3.broadPeak.bigwig")
bed2bigwig("E096-H3K4me1.broadPeak.bed","E096-H3K4me1.broadPeak.bigwig")
bed2bigwig("E096-H3K4me3.broadPeak.bed","E096-H3K4me3.broadPeak.bigwig")
bed2bigwig("E096-H3K9me3.broadPeak.bed","E096-H3K9me3.broadPeak.bigwig")

setwd("~/roadmap/liver/bed/")
#肝癌 E066
##将注释文件换成二进制
bed2bigwig("E066-H3K27ac.broadPeak.bed","E066-H3K27ac.broadPeak.bigwig")
bed2bigwig("E066-H3K27me3.broadPeak.bed","E066-H3K27me3.broadPeak.bigwig")
bed2bigwig("E066-H3K36me3.broadPeak.bed","E066-H3K36me3.broadPeak.bigwig")
bed2bigwig("E066-H3K4me1.broadPeak.bed","E066-H3K4me1.broadPeak.bigwig")
bed2bigwig("E066-H3K4me3.broadPeak.bed","E066-H3K4me3.broadPeak.bigwig")
bed2bigwig("E066-H3K9ac.broadPeak.bed","E066-H3K9ac.broadPeak.bigwig")
bed2bigwig("E066-H3K9me3.broadPeak.bed","E066-H3K9me3.broadPeak.bigwig")

setwd("~/roadmap/Pancreas/bed/")
#胰腺癌 E098
bed2bigwig("E098-DNase.hotspot.broad.bed","E098-DNase.hotspot.broad.bigwig")
bed2bigwig("E098-H3K27ac.broadPeak.bed","E098-H3K27ac.broadPeak.bigwig")
bed2bigwig("E098-H3K27me3.broadPeak.bed","E098-H3K27me3.broadPeak.bigwig")
bed2bigwig("E098-H3K36me3.broadPeak.bed","E098-H3K36me3.broadPeak.bigwig")
bed2bigwig("E098-H3K4me1.broadPeak.bed","E098-H3K4me1.broadPeak.bigwig")
bed2bigwig("E098-H3K4me3.broadPeak.bed","E098-H3K4me3.broadPeak.bigwig")
bed2bigwig("E098-H3K9me3.broadPeak.bed","E098-H3K9me3.broadPeak.bigwig")

setwd("~/roadmap/Kidney/bed/")
#肾癌 E086
bed2bigwig("E086-DNase.hotspot.broad.bed","E086-DNase.hotspot.broad.bigwig")
bed2bigwig("E086-H3K27me3.broadPeak.bed","E086-H3K27me3.broadPeak.bigwig")
bed2bigwig("E086-H3K36me3.broadPeak.bed","E086-H3K36me3.broadPeak.bigwig")
bed2bigwig("E086-H3K4me1.broadPeak.bed","E086-H3K4me1.broadPeak.bigwig")
bed2bigwig("E086-H3K4me3.broadPeak.bed","E086-H3K4me3.broadPeak.bigwig")
bed2bigwig("E086-H3K9ac.broadPeak.bed","E086-H3K9ac.broadPeak.bigwig")
bed2bigwig("E086-H3K9me3.broadPeak.bed","E086-H3K9me3.broadPeak.bigwig")

setwd("~/roadmap/Esophagus/bed/")
#食道癌 E079
bed2bigwig("E079-H3K27ac.broadPeak.bed","E079-H3K27ac.broadPeak.bigwig")
bed2bigwig("E079-H3K27me3.broadPeak.bed","E079-H3K27me3.broadPeak.bigwig")
bed2bigwig("E079-H3K36me3.broadPeak.bed","E079-H3K36me3.broadPeak.bigwig")
bed2bigwig("E079-H3K4me1.broadPeak.bed","E079-H3K4me1.broadPeak.bigwig")
bed2bigwig("E079-H3K4me3.broadPeak.bed","E079-H3K4me3.broadPeak.bigwig")
bed2bigwig("E079-H3K9me3.broadPeak.bed","E079-H3K9me3.broadPeak.bigwig")

setwd("~/roadmap/blood/bed/")
#血癌 E062
bed2bigwig("E062-H3K27ac.broadPeak.bed","E062-H3K27ac.broadPeak.bigwig")
bed2bigwig("E062-H3K27me3.broadPeak.bed","E062-H3K27me3.broadPeak.bigwig")
bed2bigwig("E062-H3K36me3.broadPeak.bed","E062-H3K36me3.broadPeak.bigwig")
bed2bigwig("E062-H3K4me1.broadPeak.bed","E062-H3K4me1.broadPeak.bigwig")
bed2bigwig("E062-H3K4me3.broadPeak.bed","E062-H3K4me3.broadPeak.bigwig")
bed2bigwig("E062-H3K9ac.broadPeak.bed","E062-H3K9ac.broadPeak.bigwig")
bed2bigwig("E062-H3K9me3.broadPeak.bed","E062-H3K9me3.broadPeak.bigwig")

setwd("~/roadmap/breast/bed/")
#乳腺癌 E028
bed2bigwig("E028-DNase.hotspot.broad.bed","E028-DNase.hotspot.broad.bigwig")
bed2bigwig("E028-H3K27me3.broadPeak.bed","E028-H3K27me3.broadPeak.bigwig")
bed2bigwig("E028-H3K36me3.broadPeak.bed","E028-H3K36me3.broadPeak.bigwig")
bed2bigwig("E028-H3K4me1.broadPeak.bed","E028-H3K4me1.broadPeak.bigwig")
bed2bigwig("E028-H3K4me3.broadPeak.bed","E028-H3K4me3.broadPeak.bigwig")
bed2bigwig("E028-H3K9me3.broadPeak.bed","E028-H3K9me3.broadPeak.bigwig")


#黑色素瘤 E059 E061
setwd("~/roadmap/E059_Melanoma/bed/")
bed2bigwig("E059-DNase.hotspot.broad.bed","E059-DNase.hotspot.broad.bigwig")
bed2bigwig("E059-H3K27ac.broadPeak.bed","E059-H3K27ac.broadPeak.bigwig")
bed2bigwig("E059-H3K27me3.broadPeak.bed","E059-H3K27me3.broadPeak.bigwig")
bed2bigwig("E059-H3K36me3.broadPeak.bed","E059-H3K36me3.broadPeak.bigwig")
bed2bigwig("E059-H3K4me1.broadPeak.bed","E059-H3K4me1.broadPeak.bigwig")
bed2bigwig("E059-H3K4me3.broadPeak.bed","E059-H3K4me3.broadPeak.bigwig")
bed2bigwig("E059-H3K9me3.broadPeak.bed","E059-H3K9me3.broadPeak.bigwig")

setwd("~/roadmap/E061_Melanoma/bed/")
bed2bigwig("E061-H3K27ac.broadPeak.bed","E061-H3K27ac.broadPeak.bigwig")
bed2bigwig("E061-H3K27me3.broadPeak.bed","E061-H3K27me3.broadPeak.bigwig")
bed2bigwig("E061-H3K36me3.broadPeak.bed","E061-H3K36me3.broadPeak.bigwig")
bed2bigwig("E061-H3K4me1.broadPeak.bed","E061-H3K4me1.broadPeak.bigwig")
bed2bigwig("E061-H3K4me3.broadPeak.bed","E061-H3K4me3.broadPeak.bigwig")
bed2bigwig("E061-H3K9me3.broadPeak.bed","E061-H3K9me3.broadPeak.bigwig")

#卵巢癌 E097
setwd("~/roadmap/Ovary/bed/")
bed2bigwig("E097-DNase.hotspot.broad.bed","E097-DNase.hotspot.broad.bigwig")
bed2bigwig("E097-H3K27ac.broadPeak.bed","E097-H3K27ac.broadPeak.bigwig")
bed2bigwig("E097-H3K27me3.broadPeak.bed","E097-H3K27me3.broadPeak.bigwig")
bed2bigwig("E097-H3K36me3.broadPeak.bed","E097-H3K36me3.broadPeak.bigwig")
bed2bigwig("E097-H3K4me1.broadPeak.bed","E097-H3K4me1.broadPeak.bigwig")
bed2bigwig("E097-H3K4me3.broadPeak.bed","E097-H3K4me3.broadPeak.bigwig")
bed2bigwig("E097-H3K9me3.broadPeak.bed","E097-H3K9me3.broadPeak.bigwig")



setwd("~/New_Analysis/")

###计算MB范围的突变数
get_reg_mution_number2 = function(orig_file, mutation_file, region_file, result_file){
  require(data.table)
  
  options(scipen = 30)
  mut <- fread(orig_file, data.table = F)
  pos_mut <- mut[mut$mut == "single base substitution",]
  mut_bed <- pos_mut[,3:5]
  mut_bed$start <- mut_bed$start - 1L
  write.table(mut_bed, mutation_file, sep = "\t", row.names = F, quote = F, col.names = F)
  
  
  rm(mut, mut_bed); gc()
  # 获取各个要计算的基因组区间
  region.coding = "merged_sorted_cds_region.bed"
  region.noncoding = "region_noncoding.bed"
  region.promoter = "merged_sorted_promoter_region.bed"
  region.nonpromoter = "region_non_promoter.bed"
  
  mutation_file.coding = file.path("~/New_Analysis/epi/mutation/", paste0("mut_coding_", basename(mutation_file)))
  mutation_file.noncoding = file.path("~/New_Analysis/epi/mutation/", paste0("mut_noncoding_", basename(mutation_file)))
  mutation_file.promoter = file.path("~/New_Analysis/epi/mutation/", paste0("mut_promoter_", basename(mutation_file)))
  mutation_file.nonpromoter = file.path("~/New_Analysis/epi/mutation/", paste0("mut_nonpromoter_", basename(mutation_file)))
  
  get_reg_mution(mutation_file, region.coding, mutation_file.coding)
  get_reg_mution(mutation_file, region.noncoding, mutation_file.noncoding)
  get_reg_mution(mutation_file, region.promoter, mutation_file.promoter)
  get_reg_mution(mutation_file, region.nonpromoter, mutation_file.nonpromoter)
  
  result_file.coding = file.path("~/New_Analysis/epi/mutation/", paste0("coding_", basename(result_file)))
  result_file.noncoding = file.path("~/New_Analysis/epi/mutation/", paste0("noncoding_", basename(result_file)))
  result_file.promoter = file.path("~/New_Analysis/epi/mutation/", paste0("promoter_", basename(result_file)))
  result_file.nonpromoter = file.path("~/New_Analysis/epi/mutation/", paste0("nonpromoter_", basename(result_file)))
  
  # 获取不同区间的1M区间突变
  get_reg_mution_number(mutation_file.coding, region_file, result_file.coding)
  get_reg_mution_number(mutation_file.noncoding, region_file, result_file.noncoding)
  get_reg_mution_number(mutation_file.promoter, region_file, result_file.promoter)
  get_reg_mution_number(mutation_file.nonpromoter, region_file, result_file.nonpromoter)
  
}

#sub("^([^.]*).*", "\\1", 'filename.extension')  remove file extension
# get_reg_mution_number2("epi/blood_cancer_info.bed", "epi/mutation/blood_cancer.bed", 
#                        "sorted_hg19_1M.bed", "1M_blood_cancer.tsv")

# 批量计算不同癌症类型和区间的突变数
for (f in dir("epi", pattern = "_info.bed", full.names = T)){
  mut_file = file.path("epi/mutation", sub("_info", "", basename(f)))
  reg_file = "sorted_hg19_1M.bed"
  result_file = paste0("1M_", sub("_info", "", basename(f)))
  
  get_reg_mution_number2(f, mut_file, reg_file, result_file)
}


###计算每个MB区间的表观遗传值
dir.create("bigwig/output")

for (f in list.files("bigwig", pattern = "^E", full.names = TRUE, include.dirs = FALSE)) {
  out = file.path("bigwig/output/", sub("([^\\.]*)\\..*", "\\1", basename(f)))
  cmd <- paste("~/bigWigAverageOverBed", f, "sorted_hg19_1M.bed", out, sep = " ")
  system(cmd)
}


rm(list = ls()); gc()

#求两列数据的平均值
cal_mean <- function(file1, file2, result_file)
{
  file1_1 <- fread(file1, data.table = F)
  file2_2 <- fread(file2, data.table = F)
  file.copy(from = file1, to = "mela")
  file.copy(from = file2, to = "mela")
  file.remove(file1)
  file.remove(file2)
  merge_result <- merge(file1_1, file2_2, by.x = "V1", by.y = "V1")
  mt <- merge_result[,c(5,10)]
  mt <- as.matrix(mt)
  mean_value <- apply(mt, 1, mean)
  merge_result$mean <- mean_value
  df <- merge_result[,c(1:4,12,6)]
  write.table(df, result_file, sep = "\t", row.names = F, col.names = F, quote = F)
}


setwd("~/New_Analysis/bigwig/output/")
dir.create("mela") # back up files not used
cal_mean("E059-H3K36me3","E061-H3K36me3","E059_E061-H3K36me3")
cal_mean("E059-H3K9me3","E061-H3K9me3","E059_E061-H3K9me3")
cal_mean("E059-H3K27ac","E061-H3K27ac","E059_E061-H3K27ac")
cal_mean("E059-H3K4me1","E061-H3K4me1","E059_E061-H3K4me1")  
cal_mean("E059-H3K27me3","E061-H3K27me3","E059_E061-H3K27me3")  
cal_mean("E059-H3K4me3","E061-H3K4me3","E059_E061-H3K4me3")      

# mapping of cancer and code
# lung E096
# liver E066
# ovary E097
# pancreas E098
# kidney E086
# esophagus E079
# blood E062
# breast E028
# melanoma E059 E061

mut_files = list.files("~/New_Analysis/epi/mutation", pattern = "1M", full.names = TRUE)
# add coverage ratio
region.coding = "merged_sorted_cds_region.bed"
region.noncoding = "region_noncoding.bed"
region.promoter = "merged_sorted_promoter_region.bed"
region.nonpromoter = "region_non_promoter.bed"

add_coverage(region.coding, mut_coding, "1M_mut_coding_add_cov.tsv")
sapply(grep("/coding", mut_files, value = TRUE), function(x) {
  add_coverage(region.coding, x, file.path(dirname(x), paste0("Cov_", basename(x))))
})

sapply(grep("/noncoding", mut_files, value = TRUE), function(x) {
  add_coverage(region.noncoding, x, file.path(dirname(x), paste0("Cov_", basename(x))))
})

sapply(grep("/promoter", mut_files, value = TRUE), function(x) {
  add_coverage(region.promoter, x, file.path(dirname(x), paste0("Cov_", basename(x))))
})

sapply(grep("/nonpromoter", mut_files, value = TRUE), function(x) {
  add_coverage(region.nonpromoter, x, file.path(dirname(x), paste0("Cov_", basename(x))))
})

mut_files = list.files("~/New_Analysis/epi/mutation", pattern = "Cov", full.names = TRUE)
annot_files = list.files("~/New_Analysis/bigwig/output/", pattern = "^E", full.names = TRUE)
mappings = c("blood", "breast", "esophagus", "kidney", "liver", "lung", "ovary", "pancreas", "mela")
names(mappings) = c("E062", "E028", "E079", "E086", "E066", "E096", "E097", "E098", "E059")


obtain_cor2 = function(mappings, mut_files, annot_files, types = c("coding", "noncoding", "promoter", "nonpromoter")) {
  require(data.table)
  options(scipen = 30)
  out = data.table()
  for (i in seq_along(mappings)) {
    message("Processing tumor type - ", mappings[i], " - ", names(mappings)[i])
    message(">> Searching mutation files...")
    type_mut = grep(mappings[i], mut_files, value = TRUE)
    message(">> OK")
    message(">> Searching annotation files...")
    type_anno = grep(names(mappings)[i], annot_files, value = TRUE)
    message(">> OK")
    for (type in types) {
      message(">> Extracting region - ", type)
      region_mut = grep(paste0("_",type), type_mut, value = TRUE) # precisely match
      message(">> OK")
      message(">> Calculating correaltion...")
      
      mut = fread(region_mut)
      colnames(mut)[6] = "cov"
      for(f in type_anno) {
        anno_type = unlist(strsplit(basename(f), split = "-"))[2]
        
        anno = fread(f)
        
        result = merge(mut, anno, by.x = "V4", by.y = "V1")
        
        # calculate correlation
        result = result[(V5.x != 0) & (V5.x != ".") & (V5.y != 0) & (V5.y != ".") ]
        result[, V5.x:=as.numeric(V5.x)]
        result[, V5.y:=as.numeric(V5.y)]
        tmp = data.table(
          tumor_type = mappings[i],
          region = type, 
          anno_type = anno_type, 
          coeff = cor(result$V5.x / result$cov, result$V5.y),
          p_val = cor.test(result$V5.x / result$cov, result$V5.y)$p.value)
        
        out = rbind(out, tmp)
      }
      message(">> OK")
    }
  }
  message("Done")
  out
}


annot_cor = obtain_cor2(mappings, mut_files, annot_files)

setwd("~/New_Analysis/")
# save
annot_cor2 = data.table::copy(annot_cor)
setnames(annot_cor2,
         colnames(annot_cor2),
         c("Tumor type" ,"Genomic region", 
           "Epigenetic type",
           "Correlation coefficient",
           "P value"))
annot_cor2[, `Genomic region` := 
                  ifelse(`Genomic region` == "nonpromoter",
                         "Noncoding-promoter", `Genomic region`)]
write.csv(as.data.frame(annot_cor2), file = "表观相关性结果2-wsx.csv", quote = FALSE, row.names = FALSE)

# plot
plot_df = list()
plot_df[["noncoding"]] = annot_cor["noncoding", on = "region"][, c("region", "p_val"):=NULL]
plot_df[["coding"]] = annot_cor["coding", on = "region"][, c("region", "p_val"):=NULL]
plot_df[["nonpromoter"]] = annot_cor["nonpromoter", on = "region"][, c("region", "p_val"):=NULL]
plot_df[["promoter"]] = annot_cor["promoter", on = "region"][, c("region", "p_val"):=NULL]

plot_heatmap = function(df){
  df_wide = tidyr::spread(data = df, tumor_type, coeff) %>% 
    tibble::column_to_rownames("anno_type")
  rownames(df_wide)[rownames(df_wide) == "DNase"] = "DNase.hotspot"
  colnames(df_wide)[colnames(df_wide) == "mela"] = "melanoma"
  # row_order = c("H3K4me3", "DNase.hotspot", "H3K9ac", 
  #               "H3K4me1", "H3K27ac", "H3K36me3",
  #               "H3K27me3", "H3K9me3")
  # col_order = c("esophagus", "lung", "kidney",
  #               "blood", "breast", "liver", 
  #               "melanoma", "ovary", "pancreas")
  pheatmap::pheatmap(df_wide, cluster_rows = F, cluster_cols = F)
}

plot_heatmap(plot_df[["noncoding"]])
plot_heatmap(plot_df[["coding"]])
plot_heatmap(plot_df[["promoter"]])
plot_heatmap(plot_df[["nonpromoter"]])

# output pdf set width 5.5, height 5 
