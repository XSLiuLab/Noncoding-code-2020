###按不同组织将肿瘤分类
mut <- fread("预处理所有突变.tsv")
donor <- fread("donor.tsv")
index <- donor[,1:2]

merge_result <- merge(mut, index, by.x = "donor", by.y = "icgc_donor_id")


#将不同的肿瘤类型输出
#肝癌
liver <- c("LIAD-FR", "LICA-CN", "LICA-FR", "LINC-JP", "LIRI-JP")
liver_cancer <- merge_result[merge_result$project_code %in% liver,]
write.table(liver_cancer, "liver_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(liver_cancer)

#胃癌
gastric <- c("STAD-US")
gastric_cancer <- merge_result[merge_result$project_code %in% gastric,]
write.table(gastric_cancer, "gastric_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(gastric_cancer)

#肾癌
kidney <- c("RECA-EU")
kidney_cancer <- merge_result[merge_result$project_code %in% kidney,]
write.table(kidney_cancer, "kidney_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(kidney_cancer)

#肺癌
lung <- c("LUSC-CN","LUSC-KR")
lung_cancer <- merge_result[merge_result$project_code %in% lung,]
write.table(lung_cancer, "lung_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(lung_cancer)

#乳腺癌
breast <- c("BRCA-EU","BRCA-FR","BRCA-US")
breast_cancer <- merge_result[merge_result$project_code %in% breast,]
write.table(breast_cancer, "breast_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(breast_cancer)

#胰腺癌
pancreas <- c("PACA-AU","PACA-CA","PAEN-AU","PAEN-IT")
pancreas_cancer <- merge_result[merge_result$project_code %in% pancreas,]
write.table(pancreas_cancer, "pancreas_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(pancreas_cancer)

#卵巢癌
ovary <- c("OV-AU")
ovary_cancer <- merge_result[merge_result$project_code %in% ovary,]
write.table(ovary_cancer, "ovary_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(ovary_cancer)

#黑色素瘤
mela <- c("MELA-AU","SKCA-BR","SKCM-US")
mela_cancer <- merge_result[merge_result$project_code %in% mela,]
write.table(mela_cancer, "mela_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(mela_cancer)

#食道癌
esophagus <- c("ESAD-UK","ESCA-CN")
esophagus_cancer <- merge_result[merge_result$project_code %in% esophagus,]
write.table(esophagus_cancer, "esophagus_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(esophagus_cancer)

#血癌
blood <- c("ALL-US","CLLE-ES","MALY-DE","NKTL-SG")
blood_cancer <- merge_result[merge_result$project_code %in% blood,]
write.table(blood_cancer, "blood_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(blood_cancer)

#结肠癌
colon <- c("COAD-US","COCA-CN")
colon_cancer <- merge_result[merge_result$project_code %in% colon,]
write.table(colon_cancer, "colon_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(colon_cancer)

#胶质瘤
brain <- c("LGG-US")
brain_cancer <- merge_result[merge_result$project_code %in% brain,]
write.table(brain_cancer, "brain_cancer_info.bed", sep = "\t", row.names = F, quote = F)
rm(brain_cancer)









##将表观遗传bed文件转化为bigwig文件
bed2bigwig <- function(bedfile, result_file)
{
		temp1 <- paste("sorted_", bedfile, sep = "")
		cmd1 <- paste("sort -k 1,1 -k 2,2n ",bedfile," > ",temp1, sep = "")
		system(cmd1)
		cmd2 <- paste("~/bedGraphToBigWig ",temp1, "~/hg19_chr_size.bed ",result_file)
		system(cmd2)
		file.remove(temp1)
}


#删除chrM染色体
rm_chrM <- function(file)
{
	df <- fread(file, data.table = F)
	df <- df[df$V1 != "chrM",]
	write.table(df, file, sep = "\t", row.names = F, quote = F, col.names = F)
}

#将DNase文件转换为bed 格式
toBed <- function(file)
{
	df <- fread(file, data.table = F)
	df <- df[,c(1:3,5)]
	write.table(df, file, sep = "\t", row.names = F, quote = F, col.names = F)
}



#肺癌 E096
##将注释文件换成二进制
bed2bigwig("E096-H3K27ac.broadPeak.bed","E096-H3K27ac.broadPeak.bigwig")
bed2bigwig("E096-H3K27me3.broadPeak.bed","E096-H3K27me3.broadPeak.bigwig")
bed2bigwig("E096-H3K36me3.broadPeak.bed","E096-H3K36me3.broadPeak.bigwig")
bed2bigwig("E096-H3K4me1.broadPeak.bed","E096-H3K4me1.broadPeak.bigwig")
bed2bigwig("E096-H3K4me3.broadPeak.bed","E096-H3K4me3.broadPeak.bigwig")
bed2bigwig("E096-H3K9me3.broadPeak.bed","E096-H3K9me3.broadPeak.bigwig")


#肝癌 E066
##将注释文件换成二进制
bed2bigwig("E066-H3K27ac.broadPeak.bed","E066-H3K27ac.broadPeak.bigwig")
bed2bigwig("E066-H3K27me3.broadPeak.bed","E066-H3K27me3.broadPeak.bigwig")
bed2bigwig("E066-H3K36me3.broadPeak.bed","E066-H3K36me3.broadPeak.bigwig")
bed2bigwig("E066-H3K4me1.broadPeak.bed","E066-H3K4me1.broadPeak.bigwig")
bed2bigwig("E066-H3K4me3.broadPeak.bed","E066-H3K4me3.broadPeak.bigwig")
bed2bigwig("E066-H3K9ac.broadPeak.bed","E066-H3K9ac.broadPeak.bigwig")
bed2bigwig("E066-H3K9me3.broadPeak.bed","E066-H3K9me3.broadPeak.bigwig")


#卵巢癌 E097
##将注释文件换成二进制
bed2bigwig("E097-DNase.hotspot.broad.bed","E097-DNase.hotspot.broad.bigwig")
bed2bigwig("E097-H3K27ac.broadPeak.bed","E097-H3K27ac.broadPeak.bigwig")
bed2bigwig("E097-H3K27me3.broadPeak.bed","E097-H3K27me3.broadPeak.bigwig")
bed2bigwig("E097-H3K36me3.broadPeak.bed","E097-H3K36me3.broadPeak.bigwig")
bed2bigwig("E097-H3K4me1.broadPeak.bed","E097-H3K4me1.broadPeak.bigwig")
bed2bigwig("E097-H3K4me3.broadPeak.bed","E097-H3K4me3.broadPeak.bigwig")
bed2bigwig("E097-H3K9me3.broadPeak.bed","E097-H3K9me3.broadPeak.bigwig")

   
#胃癌 E094
bed2bigwig("E094-DNase.hotspot.broad.bed","E094-DNase.hotspot.broad.bigwig")
bed2bigwig("E094-H3K27ac.broadPeak.bed","E094-H3K27ac.broadPeak.bigwig")
bed2bigwig("E094-H3K27me3.broadPeak.bed","E094-H3K27me3.broadPeak.bigwig")
bed2bigwig("E094-H3K36me3.broadPeak.bed","E094-H3K36me3.broadPeak.bigwig")
bed2bigwig("E094-H3K4me1.broadPeak.bed","E094-H3K4me1.broadPeak.bigwig")
bed2bigwig("E094-H3K4me3.broadPeak.bed","E094-H3K4me3.broadPeak.bigwig")
bed2bigwig("E094-H3K9me3.broadPeak.bed","E094-H3K9me3.broadPeak.bigwig")


#胰腺癌 E098
bed2bigwig("E098-DNase.hotspot.broad.bed","E098-DNase.hotspot.broad.bigwig")
bed2bigwig("E098-H3K27ac.broadPeak.bed","E098-H3K27ac.broadPeak.bigwig")
bed2bigwig("E098-H3K27me3.broadPeak.bed","E098-H3K27me3.broadPeak.bigwig")
bed2bigwig("E098-H3K36me3.broadPeak.bed","E098-H3K36me3.broadPeak.bigwig")
bed2bigwig("E098-H3K4me1.broadPeak.bed","E098-H3K4me1.broadPeak.bigwig")
bed2bigwig("E098-H3K4me3.broadPeak.bed","E098-H3K4me3.broadPeak.bigwig")
bed2bigwig("E098-H3K9me3.broadPeak.bed","E098-H3K9me3.broadPeak.bigwig")


#肾癌 E086
bed2bigwig("E086-DNase.hotspot.broad.bed","E086-DNase.hotspot.broad.bigwig")
bed2bigwig("E086-H3K27me3.broadPeak.bed","E086-H3K27me3.broadPeak.bigwig")
bed2bigwig("E086-H3K36me3.broadPeak.bed","E086-H3K36me3.broadPeak.bigwig")
bed2bigwig("E086-H3K4me1.broadPeak.bed","E086-H3K4me1.broadPeak.bigwig")
bed2bigwig("E086-H3K4me3.broadPeak.bed","E086-H3K4me3.broadPeak.bigwig")
bed2bigwig("E086-H3K9ac.broadPeak.bed","E086-H3K9ac.broadPeak.bigwig")
bed2bigwig("E086-H3K9me3.broadPeak.bed","E086-H3K9me3.broadPeak.bigwig")


#食道癌 E079
bed2bigwig("E079-H3K27ac.broadPeak.bed","E079-H3K27ac.broadPeak.bigwig")
bed2bigwig("E079-H3K27me3.broadPeak.bed","E079-H3K27me3.broadPeak.bigwig")
bed2bigwig("E079-H3K36me3.broadPeak.bed","E079-H3K36me3.broadPeak.bigwig")
bed2bigwig("E079-H3K4me1.broadPeak.bed","E079-H3K4me1.broadPeak.bigwig")
bed2bigwig("E079-H3K4me3.broadPeak.bed","E079-H3K4me3.broadPeak.bigwig")
bed2bigwig("E079-H3K9me3.broadPeak.bed","E079-H3K9me3.broadPeak.bigwig")


#血癌 E062
bed2bigwig("E062-H3K27ac.broadPeak.bed","E062-H3K27ac.broadPeak.bigwig")
bed2bigwig("E062-H3K27me3.broadPeak.bed","E062-H3K27me3.broadPeak.bigwig")
bed2bigwig("E062-H3K36me3.broadPeak.bed","E062-H3K36me3.broadPeak.bigwig")
bed2bigwig("E062-H3K4me1.broadPeak.bed","E062-H3K4me1.broadPeak.bigwig")
bed2bigwig("E062-H3K4me3.broadPeak.bed","E062-H3K4me3.broadPeak.bigwig")
bed2bigwig("E062-H3K9ac.broadPeak.bed","E062-H3K9ac.broadPeak.bigwig")
bed2bigwig("E062-H3K9me3.broadPeak.bed","E062-H3K9me3.broadPeak.bigwig")


#乳腺癌 E028
bed2bigwig("E028-DNase.hotspot.broad.bed","E028-DNase.hotspot.broad.bigwig")
bed2bigwig("E028-H3K27me3.broadPeak.bed","E028-H3K27me3.broadPeak.bigwig")
bed2bigwig("E028-H3K36me3.broadPeak.bed","E028-H3K36me3.broadPeak.bigwig")
bed2bigwig("E028-H3K4me1.broadPeak.bed","E028-H3K4me1.broadPeak.bigwig")
bed2bigwig("E028-H3K4me3.broadPeak.bed","E028-H3K4me3.broadPeak.bigwig")
bed2bigwig("E028-H3K9me3.broadPeak.bed","E028-H3K9me3.broadPeak.bigwig")


#胶质瘤 E053
bed2bigwig("E053-H3K27me3.broadPeak.bed","E053-H3K27me3.broadPeak.bigwig")
bed2bigwig("E053-H3K36me3.broadPeak.bed","E053-H3K36me3.broadPeak.bigwig")
bed2bigwig("E053-H3K4me1.broadPeak.bed","E053-H3K4me1.broadPeak.bigwig")
bed2bigwig("E053-H3K4me3.broadPeak.bed","E053-H3K4me3.broadPeak.bigwig")
bed2bigwig("E053-H3K9me3.broadPeak.bed","E053-H3K9me3.broadPeak.bigwig")


#结肠癌 E106
bed2bigwig("E106-H3K27ac.broadPeak.bed","E106-H3K27ac.broadPeak.bigwig")
bed2bigwig("E106-H3K27me3.broadPeak.bed","E106-H3K27me3.broadPeak.bigwig")
bed2bigwig("E106-H3K36me3.broadPeak.bed","E106-H3K36me3.broadPeak.bigwig")
bed2bigwig("E106-H3K4me3.broadPeak.bed","E106-H3K4me3.broadPeak.bigwig")
bed2bigwig("E106-H3K4me1.broadPeak.bed","E106-H3K4me1.broadPeak.bigwig")
bed2bigwig("E106-H3K9me3.broadPeak.bed","E106-H3K9me3.broadPeak.bigwig")


#黑色素瘤 E059 E061
bed2bigwig("E059-DNase.hotspot.broad.bed","E059-DNase.hotspot.broad.bigwig")
bed2bigwig("E059-H3K27ac.broadPeak.bed","E059-H3K27ac.broadPeak.bigwig")
bed2bigwig("E059-H3K27me3.broadPeak.bed","E059-H3K27me3.broadPeak.bigwig")
bed2bigwig("E059-H3K36me3.broadPeak.bed","E059-H3K36me3.broadPeak.bigwig")
bed2bigwig("E059-H3K4me1.broadPeak.bed","E059-H3K4me1.broadPeak.bigwig")
bed2bigwig("E059-H3K4me3.broadPeak.bed","E059-H3K4me3.broadPeak.bigwig")
bed2bigwig("E059-H3K9me3.broadPeak.bed","E059-H3K9me3.broadPeak.bigwig")


bed2bigwig("E061-H3K27ac.broadPeak.bed","E061-H3K27ac.broadPeak.bigwig")
bed2bigwig("E061-H3K27me3.broadPeak.bed","E061-H3K27me3.broadPeak.bigwig")
bed2bigwig("E061-H3K36me3.broadPeak.bed","E061-H3K36me3.broadPeak.bigwig")
bed2bigwig("E061-H3K4me1.broadPeak.bed","E061-H3K4me1.broadPeak.bigwig")
bed2bigwig("E061-H3K4me3.broadPeak.bed","E061-H3K4me3.broadPeak.bigwig")
bed2bigwig("E061-H3K9me3.broadPeak.bed","E061-H3K9me3.broadPeak.bigwig")



 
 
 
###计算MB范围的非编码突变数

#改为单碱基突变的bed格式
mut2bed <- function(file_names)
{
	mut <- fread(file_names, data.table = F)
	pos_mut <- mut[mut$type == "single base substitution",]
	mut_bed <- pos_mut[,3:5]
	mut_bed$start <- mut_bed$start - 1
	write.table(mut_bed, file_names, sep = "\t", row.names = F, quote = F, col.names = F)
}

#计算hg19 MB 范围内的突变次数
count_mb <- function(file_names, result_file)
{
	mut <- fread(file_names, data.table = F)
	mut$V4 <- 1
	write.table(mut, file_names, sep = "\t", row.names = F, col.names = F, quote = F)
	
	cmd <- paste("bedtools map -a  ~/sort_hg19_mb.bed -b ",file_names, " -c 4 > ", result_file, sep = "")
	system(cmd)
}

#删除科学计数法
rm_e <- function(file_names)
{
	options(scipen = 30)
	a <- fread(file_names, data.table = F)
	write.table(a, file_names, sep = "\t", row.names = F, col.names = F, quote = F)
}



#血癌
mut2bed("blood_cancer_info.bed")
#选取非编码区突变
rm_e("blood_cancer_info.bed")
region_mut("blood_cancer_info.bed","hg19_noncoding_region.bed","noncod_blood_mut.bed")
#计算MB非编码突变数
count_mb("noncod_blood_mut.bed","noncod_blood_MBnum.bed")

#胶质瘤（突变数目太少，不做）
mut2bed("brain_cancer_info.bed")
region_mut("brain_cancer_info.bed","hg19_noncoding_region.bed","noncod_brain_mut.bed")
count_mb("noncod_brain_mut.bed","noncod_brain_MBnum.bed")

#乳腺癌
mut2bed("breast_cancer_info.bed")
region_mut("breast_cancer_info.bed","hg19_noncoding_region.bed","noncod_breast_mut.bed")
count_mb("noncod_breast_mut.bed","noncod_breast_MBnum.bed")

#结肠癌
mut2bed("colon_cancer_info.bed")
region_mut("colon_cancer_info.bed","hg19_noncoding_region.bed","noncod_colon_mut.bed")
count_mb("noncod_colon_mut.bed","noncod_colon_MBnum.bed")

#食管癌
mut2bed("esophagus_cancer_info.bed")
region_mut("esophagus_cancer_info.bed","hg19_noncoding_region.bed","noncod_esophagus_mut.bed")
count_mb("noncod_esophagus_mut.bed","noncod_esophagus_MBnum.bed")

#肾癌
mut2bed("kidney_cancer_info.bed")
region_mut("kidney_cancer_info.bed","hg19_noncoding_region.bed","noncod_kidney_mut.bed")
count_mb("noncod_kidney_mut.bed","noncod_kidney_MBnum.bed")

#肝癌
mut2bed("liver_cancer_info.bed")
region_mut("liver_cancer_info.bed","hg19_noncoding_region.bed","noncod_liver_mut.bed")
count_mb("noncod_liver_mut.bed","noncod_liver_MBnum.bed")

#肺癌
mut2bed("lung_cancer_info.bed")
region_mut("lung_cancer_info.bed","hg19_noncoding_region.bed","noncod_lung_mut.bed")
count_mb("noncod_lung_mut.bed","noncod_lung_MBnum.bed")

#皮肤癌
mut2bed("mela_cancer_info.bed")
region_mut("mela_cancer_info.bed","hg19_noncoding_region.bed","noncod_mela_mut.bed")
count_mb("noncod_mela_mut.bed","noncod_mela_MBnum.bed")

#卵巢癌
mut2bed("ovary_cancer_info.bed")
region_mut("ovary_cancer_info.bed","hg19_noncoding_region.bed","noncod_ovary_mut.bed")
count_mb("noncod_ovary_mut.bed","noncod_ovary_MBnum.bed")

#胰腺癌
mut2bed("pancreas_cancer_info.bed")
region_mut("pancreas_cancer_info.bed","hg19_noncoding_region.bed","noncod_pancreas_mut.bed")
count_mb("noncod_pancreas_mut.bed","noncod_pancreas_MBnum.bed")

#胃癌（突变数目太少，不做）
mut2bed("gastric_cancer_info.bed")
region_mut("gastric_cancer_info.bed","hg19_noncoding_region.bed","noncod_gastric_mut.bed")
count_mb("noncod_gastric_mut.bed","noncod_gastric_MBnum.bed")




###计算每个MB区间的表观遗传值
library(stringr)

#对文件夹中的每个文件进行计算
cal_all <- function()
{
	file_names <- list.files()
	for(i in file_names)
	{
		anno_type <- str_split_fixed(i, ".broadPeak.bigwig", n = 2)[,1] 
		cmd <- paste("~/bigWigAverageOverBed ",i," ~/sort_hg19_mb.bed ",anno_type, sep = "")
		system(cmd)
	}
}


#血癌
setwd("/home/zhangjing/paper/epi_cor/blood/E062/bigwig")
cal_all()

#乳腺癌
setwd("/home/zhangjing/paper/epi_cor/breast/E028/bigwig")
cal_all()

#结直肠癌
setwd("/home/zhangjing/paper/epi_cor/colon/E106/bigwig")
cal_all()

#食道癌
setwd("/home/zhangjing/paper/epi_cor/esophagus/E079/bigwig")
cal_all()

#肾癌
setwd("/home/zhangjing/paper/epi_cor/kidney/E086/bigwig")
cal_all()

#肝癌
setwd("/home/zhangjing/paper/epi_cor/liver/E066/bigwig")
cal_all()

#肺癌
setwd("/home/zhangjing/paper/epi_cor/lung/E096/bigwig")
cal_all()

#卵巢癌
setwd("/home/zhangjing/paper/epi_cor/ovary/E097/bigwig")
cal_all()

#胰腺癌
setwd("/home/zhangjing/paper/epi_cor/pancreas/E098/bigwig")
cal_all()

#皮肤癌一
setwd("/home/zhangjing/paper/epi_cor/mela/E059/bigwig")
cal_all()
#皮肤癌二
setwd("/home/zhangjing/paper/epi_cor/mela/E061/bigwig")
cal_all()





###计算相关性
cal_cor <- function()
{
	file_names <- list.files()
	mut_file <- file_names[grep("MBnum.bed", file_names)]
	anno_file <- file_names[-grep("MBnum.bed", file_names)]
	
	mut <- fread(mut_file, data.table = F)
	mut <- mut[mut$V5 != ".",]
	mut$V5 <- as.numeric(mut$V5)
	
	for(i in anno_file)
	{
		anno <- fread(i, data.table = F)
		anno <- anno[anno$V5 != 0,]
		result <- merge(mut, anno, by.x = "V4", by.y = "V1")
		cor <- cor(result$V5.x, result$V5.y)
		p <- cor.test(result$V5.x, result$V5.y)$p.value
		res1 <- paste(i,"cor is ", cor, sep = "  ")
		res2 <- paste(i,"p value is ", p, sep = "  ")
		print(res1)
		print(res2)
	}
	
}

#求两列数据的平均值
cal_mean <- function(file1, file2, result_file)
{
	file1 <- fread(file1, data.table = F)
	file2 <- fread(file2, data.table = F)
	merge_result <- merge(file1, file2, by.x = "V1", by.y = "V1")
	mt <- merge_result[,c(5,10)]
	mt <- as.matrix(mt)
	mean_value <- apply(mt, 1, mean)
	merge_result$mean <- mean_value
	df <- merge_result[,c(1:4,12,6)]
	write.table(df, result_file, sep = "\t", row.names = F, col.names = F, quote = F)
}





#血癌
setwd("/home/zhangjing/paper/epi_cor/blood/cor")
cal_cor()

#乳腺癌
setwd("/home/zhangjing/paper/epi_cor/breast/cor")
cal_cor()

#结直肠癌
setwd("/home/zhangjing/paper/epi_cor/colon/cor")
cal_cor()

#食道癌
setwd("/home/zhangjing/paper/epi_cor/esophagus/cor/")
cal_cor()

#肾癌
setwd("/home/zhangjing/paper/epi_cor/kidney/cor")
cal_cor()

#肝癌
setwd("/home/zhangjing/paper/epi_cor/liver/cor")
cal_cor()

#肺癌
setwd("/home/zhangjing/paper/epi_cor/lung/cor")
cal_cor()

#卵巢癌
setwd("/home/zhangjing/paper/epi_cor/ovary/cor")
cal_cor()

#胰腺癌
setwd("/home/zhangjing/paper/epi_cor/pancreas/cor")
cal_cor()


#皮肤癌（对两个数据取平均值）
setwd("~/paper/epi_cor/mela/cor/")

cal_mean("E059-H3K36me3","E061-H3K36me3","E059_E061-H3K36me3")
cal_mean("E059-H3K9me3","E061-H3K9me3","E059_E061-H3K9me3")
cal_mean("E059-H3K27ac","E061-H3K27ac","E059_E061-H3K27ac")
cal_mean("E059-H3K4me1","E061-H3K4me1","E059_E061-H3K4me1")  
cal_mean("E059-H3K27me3","E061-H3K27me3","E059_E061-H3K27me3")  
cal_mean("E059-H3K4me3","E061-H3K4me3","E059_E061-H3K4me3")      
       
cal_cor()
   	   

 
   
  











