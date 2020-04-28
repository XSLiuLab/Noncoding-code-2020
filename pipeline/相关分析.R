###将hg19基因组分割为不同区间的片段
#human_genome_chromosome_length.tsv为1:22+x+y染色体的长度，从human.hg19.genome文件中获取


get_region_bed <- function(region_length)
{
  hg19_chr <- read.table("human_genome_chromosome_length.tsv", header = F)
  
  chr <- as.character(unique(hg19_chr$V1))
  chr_list <- vector("list", length = length(chr))
  names(chr_list)  <- chr
  
  for(i in chr)
  {
    end_pos <- hg19_chr[hg19_chr$V1 == i,]$V2
    df_end <- seq(0, end_pos, region_length)[-1]
    df_start <- df_end - region_length
    df_chr <- as.character(hg19_chr[hg19_chr$V1 == i,]$V1)
    df <- data.frame(df_chr, df_start, df_end)
    chr_list[[i]] <- df
  }
  
  result <- do.call("rbind", chr_list)
  result$name <- 1:dim(result)[1]
  write.table(result, "hg19_mb(500kb/100kb/10kb/1kb).bed", sep = "\t", row.names = F, quote = F, col.names = F)
}  



######################################计算不同区间上的突变数目######################################
options(scipen = 30)



###获取基因组上所有的非编码区
system("sort -k 1,1 -k 2,2n cds_region.bed > sorted_cds_region.bed")   #前面分析中得到cds_region.bed
system("sort -k 1,1 -k 2,2n human.hg19.genome > sorted_human.hg19.genome")  #bedtools包中附带文件
system("bedtools complement -i sorted_cds_region.bed -g sorted_human.hg19.genome > hg19_noncoding_region.bed")


###获取非编码区上的所有突变
mut <- fread("预处理所有突变.tsv", data.table = F)
pos_mut <- mut[mut$type == "single base substitution",]
###以bed格式输出所有单碱基突变 
mut_bed <- pos_mut[,3:5]
mut_bed$start <- mut_bed$start - 1
write.table(mut_bed, "single_mut.bed", sep = "\t", row.names = F, quote = F, col.names = F)

region_mut("single_mut.bed","hg19_noncoding_region.bed","all_noncoding_mut")
file.remove("cds_region.bed","sorted_cds_region.bed","sorted_human.hg19.genome","hg19_noncoding_region.bed","single_mut.bed")



###将每个突变后面添加数字1
mut <- fread("all_noncoding_mut", data.table = F)
mut$V4 <- 1
write.table(mut,"all_noncoding_mut", sep = "\t", row.names = F, col.names = F, quote = F)

###计算每MB区间内的突变次数
system("bedtools map -a sort_hg19_mb.bed -b all_noncoding_mut -c 4 > hg19_noncod_mut_numberMB.tsv")


###删除缺失值
mut_num <- read.table("hg19_noncod_mut_numberMB.tsv", stringsAsFactors = F)
mut_num <- mut_num[mut_num$V5 != ".",]
mut_num$V5 <- as.numeric(mut_num$V5)



#####将每个1mb的突变数目和相对应的注释特征匹配
#重组率
rec_rata <- read.table("rec_rate.tsv", header = T)
result <- merge(rec_rata, mut_num, by.x = "V4", by.y = "V4")


#拼接性
map <- read.table("mean_mb_mappability.tsv")
result <- merge(map, mut_num, by.x = "V1", by.y = "V4")

#保守性
con <- read.table("conservation_Mb.bed")
result <- merge(con, mut_num, by.x = "V1", by.y = "V4")

#复制时间
reptime <- read.table("reptime_Mb.bed", stringsAsFactors = F)
reptime <- reptime[reptime$V5 != ".",]
reptime$V5 <- as.numeric(reptime$V5)
result <- merge(reptime, mut_num, by.x = "V4", by.y = "V4")

#TFBS
tfbs <- read.table("tfbs_mb.bed")
result <- merge(tfbs, mut_num, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]

#GC含量
gc <- read.table("hg19_gcMb.bed", stringsAsFactors = F)
gc <- gc[gc$V5 != ".",]
gc$V5 <- as.numeric(gc$V5)
result <- merge(gc, mut_num, by.x = "V4", by.y = "V4")

#CpG岛
cpg <- read.table("CpG_percent_Mb.bed")
result <- merge(cpg, mut_num, by.x = "V4", by.y = "V4")

#聚合酶II
pol <- read.table("polII_mb.bed")
result <- merge(pol, mut_num, by.x = "V1", by.y = "V4")





#####对每个特征画散点图
###区间非编码突变频数
mut_num <- read.table("hg19_noncod_mut_numberMB.tsv", stringsAsFactors = F)
mut_num <- mut_num[mut_num$V5 != ".",]
mut_num$V5 <- as.numeric(mut_num$V5)



#拼接性
map <- read.table("mean_mb_mappability.tsv")
result <- merge(map, mut_num, by.x = "V1", by.y = "V4")
smoothScatter(result$V5.x, result$V5.y)


#保守性
con <- read.table("conservation_Mb.bed")
result <- merge(con, mut_num, by.x = "V1", by.y = "V4")
smoothScatter(result$V5.x, result$V5.y)


#转录因子结合位点
tfbs <- read.table("tfbs_mb.bed")
result <- merge(tfbs, mut_num, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
smoothScatter(result$V5.x, result$V5.y)



#复制时间
reptime <- read.table("reptime_Mb.bed", stringsAsFactors = F)
reptime <- reptime[reptime$V5 != ".",]
reptime$V5 <- as.numeric(reptime$V5)
result <- merge(reptime, mut_num, by.x = "V4", by.y = "V4")
smoothScatter(result$V5.x, result$V5.y)



#GC含量
gc <- read.table("hg19_gcMb.bed", stringsAsFactors = F)
gc <- gc[gc$V5 != ".",]
gc$V5 <- as.numeric(gc$V5)
result <- merge(gc, mut_num, by.x = "V4", by.y = "V4")
smoothScatter(result$V5.x, result$V5.y)



#CpG岛
cpg <- read.table("CpG_percent_Mb.bed")
result <- merge(cpg, mut_num, by.x = "V4", by.y = "V4")
smoothScatter(result$V5.x, result$V5.y)



#重组率相关系数
rec_rata <- read.table("rec_rate.tsv", header = T)
result <- merge(rec_rata, mut_num, by.x = "V4", by.y = "V4")
smoothScatter(result$avg, result$V5.y)



#聚合酶II
pol <- read.table("polII_mb.bed")
result <- merge(pol, mut_num, by.x = "V1", by.y = "V4")
smoothScatter(result$V3.x, result$V5)






#########获取不同区间大小的保守性数值########
{
  ./bigWigAverageOverBed hg19.100way.phastCons.bw sort_hg19_mb.bed conservation_Mb.bed
  ./bigWigAverageOverBed hg19.100way.phastCons.bw sort_hg19_500kb.bed conservation_500kb.bed
  ./bigWigAverageOverBed hg19.100way.phastCons.bw sort_hg19_100kb.bed conservation_100kb.bed
  ./bigWigAverageOverBed hg19.100way.phastCons.bw sort_hg19_10kb.bed conservation_10kb.bed
}


##########计算每个区间内CpG岛的长度#########
{
  bedtools map -a sort_hg19_mb.bed -b  sort_cpg_length.bed  -c 4 -o sum >  CpG_length_sum_Mb.bed
  bedtools map -a sort_hg19_500kb.bed -b  sort_cpg_length.bed  -c 4 -o sum >  CpG_length_sum_500kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  sort_cpg_length.bed  -c 4 -o sum >  CpG_length_sum_100kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  sort_cpg_length.bed  -c 4 -o sum >  CpG_length_sum_10kb.bed
  
  
  files <- c("CpG_length_sum_Mb.bed","CpG_length_sum_500kb.bed","CpG_length_sum_100kb.bed","CpG_length_sum_10kb.bed")
  read_file <- vector("list", length = 4)
  names(read_file) <- files
  for(i in files) read_file[[i]] <- fread(i)
  
  for(i in files)
  {
    df <- read_file[[i]]
    df <- df[df$V5 != ".",]
    df$V5 <- as.numeric(df$V5)
    read_file[[i]] <- df
  }
  
  read_file$CpG_length_sum_Mb.bed$V5 <- read_file$CpG_length_sum_Mb.bed$V5/1000000
  read_file$CpG_length_sum_500kb.bed$V5 <- read_file$CpG_length_sum_500kb.bed$V5/500000
  read_file$CpG_length_sum_100kb.bed$V5 <- read_file$CpG_length_sum_100kb.bed$V5/100000
  read_file$CpG_length_sum_10kb.bed$V5 <- read_file$CpG_length_sum_10kb.bed$V5/10000
  
  
  write.table(read_file$CpG_length_sum_Mb.bed,"CpG_percent_Mb.bed",sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(read_file$CpG_length_sum_500kb.bed,"CpG_percent_500kb.bed",sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(read_file$CpG_length_sum_100kb.bed,"CpG_percent_100kb.bed",sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(read_file$CpG_length_sum_10kb.bed,"CpG_percent_10kb.bed",sep = "\t", row.names = F, col.names = F, quote = F)
}


###########GC含量###########################
{
  library(data.table)
  library(stringr)
  options(scipen = 30)
  
  gc <- fread("hg19.gc5Base.txt", nThread = 8, sep = ",", header = F)
  chr_index <- grep("variableStep chrom=", gc$V1)
  chr_index_n <- chr_index[1:25]
  
  gc_chr <- vector("list", length = 24)
  chr_names <- gc$V1[chr_index_n][1:24]
  names(gc_chr) <- chr_names
  for(i in 1:24)   gc_chr[[chr_names[i]]] <- gc[chr_index_n[i]:chr_index_n[i+1],]
  
  for(i in chr_names)
  {
    df <- gc_chr[[i]]
    dim_length <- dim(df)[1]
    df_value <- df$V1[2:(dim_length-1)]
    pos <- str_split_fixed(df_value,"\t",2)[,1]
    gc_value <- str_split_fixed(df_value,"\t",2)[,2]
    pos_start <- as.numeric(pos) - 1
    pos_end <- as.numeric(pos) + 4
    chr <- strsplit(i,"variableStep chrom=| span=5")[[1]][2]
    bed_df <- data.frame(chr, pos_start, pos_end, gc_value)
    bed_list[[i]] <- bed_df
    fwrite(bed_df, paste(chr,"_gc5base.bed", sep = ""), sep = "\t", row.names = F, quote = F, col.names = F  )
  }
  
  
  #对 Mb 区间计算GC平均值（开多个终端，同时运算）
  bedtools map -a sort_hg19_mb.bed -b  chr1_gc5base.bed  -c 4 -o mean >  chr1_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr2_gc5base.bed  -c 4 -o mean >  chr2_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr3_gc5base.bed  -c 4 -o mean >  chr3_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr4_gc5base.bed  -c 4 -o mean >  chr4_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr5_gc5base.bed  -c 4 -o mean >  chr5_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr6_gc5base.bed  -c 4 -o mean >  chr6_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr7_gc5base.bed  -c 4 -o mean >  chr7_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr8_gc5base.bed  -c 4 -o mean >  chr8_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr9_gc5base.bed  -c 4 -o mean >  chr9_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr10_gc5base.bed  -c 4 -o mean >  chr10_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr11_gc5base.bed  -c 4 -o mean >  chr11_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr12_gc5base.bed  -c 4 -o mean >  chr12_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr13_gc5base.bed  -c 4 -o mean >  chr13_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr14_gc5base.bed  -c 4 -o mean >  chr14_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr15_gc5base.bed  -c 4 -o mean >  chr15_gcMb.bed
  
  
  bedtools map -a sort_hg19_mb.bed -b  chr16_gc5base.bed  -c 4 -o mean >  chr16_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr17_gc5base.bed  -c 4 -o mean >  chr17_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr18_gc5base.bed  -c 4 -o mean >  chr18_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr19_gc5base.bed  -c 4 -o mean >  chr19_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr20_gc5base.bed  -c 4 -o mean >  chr20_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chr21_gc5base.bed  -c 4 -o mean >  chr21_gcMb.bed
  
  bedtools map -a sort_hg19_mb.bed -b  chr22_gc5base.bed  -c 4 -o mean >  chr22_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chrX_gc5base.bed  -c 4 -o mean >  chrX_gcMb.bed
  bedtools map -a sort_hg19_mb.bed -b  chrY_gc5base.bed  -c 4 -o mean >  chrY_gcMb.bed
  
  
  #对 100kb 区间计算GC平均值
  bedtools map -a sort_hg19_100kb.bed -b  chr1_gc5base.bed  -c 4 -o mean >  chr1_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr2_gc5base.bed  -c 4 -o mean >  chr2_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr3_gc5base.bed  -c 4 -o mean >  chr3_gc100kb.bed
  
  bedtools map -a sort_hg19_100kb.bed -b  chr4_gc5base.bed  -c 4 -o mean >  chr4_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr5_gc5base.bed  -c 4 -o mean >  chr5_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr6_gc5base.bed  -c 4 -o mean >  chr6_gc100kb.bed
  
  bedtools map -a sort_hg19_100kb.bed -b  chr7_gc5base.bed  -c 4 -o mean >  chr7_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr8_gc5base.bed  -c 4 -o mean >  chr8_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr9_gc5base.bed  -c 4 -o mean >  chr9_gc100kb.bed
  
  bedtools map -a sort_hg19_100kb.bed -b  chr10_gc5base.bed  -c 4 -o mean >  chr10_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr11_gc5base.bed  -c 4 -o mean >  chr11_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr12_gc5base.bed  -c 4 -o mean >  chr12_gc100kb.bed
  
  bedtools map -a sort_hg19_100kb.bed -b  chr13_gc5base.bed  -c 4 -o mean >  chr13_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr14_gc5base.bed  -c 4 -o mean >  chr14_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr15_gc5base.bed  -c 4 -o mean >  chr15_gc100kb.bed
  
  bedtools map -a sort_hg19_100kb.bed -b  chr16_gc5base.bed  -c 4 -o mean >  chr16_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr17_gc5base.bed  -c 4 -o mean >  chr17_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr18_gc5base.bed  -c 4 -o mean >  chr18_gc100kb.bed
  
  bedtools map -a sort_hg19_100kb.bed -b  chr19_gc5base.bed  -c 4 -o mean >  chr19_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr20_gc5base.bed  -c 4 -o mean >  chr20_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chr21_gc5base.bed  -c 4 -o mean >  chr21_gc100kb.bed
  
  bedtools map -a sort_hg19_100kb.bed -b  chr22_gc5base.bed  -c 4 -o mean >  chr22_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chrX_gc5base.bed  -c 4 -o mean >  chrX_gc100kb.bed
  bedtools map -a sort_hg19_100kb.bed -b  chrY_gc5base.bed  -c 4 -o mean >  chrY_gc100kb.bed
  
  #对 10kb 区间计算GC平均值
  bedtools map -a sort_hg19_10kb.bed -b  chr1_gc5base.bed  -c 4 -o mean >  chr1_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr2_gc5base.bed  -c 4 -o mean >  chr2_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr3_gc5base.bed  -c 4 -o mean >  chr3_gc10kb.bed
  
  bedtools map -a sort_hg19_10kb.bed -b  chr4_gc5base.bed  -c 4 -o mean >  chr4_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr5_gc5base.bed  -c 4 -o mean >  chr5_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr6_gc5base.bed  -c 4 -o mean >  chr6_gc10kb.bed
  
  bedtools map -a sort_hg19_10kb.bed -b  chr7_gc5base.bed  -c 4 -o mean >  chr7_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr8_gc5base.bed  -c 4 -o mean >  chr8_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr9_gc5base.bed  -c 4 -o mean >  chr9_gc10kb.bed
  
  bedtools map -a sort_hg19_10kb.bed -b  chr10_gc5base.bed  -c 4 -o mean >  chr10_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr11_gc5base.bed  -c 4 -o mean >  chr11_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr12_gc5base.bed  -c 4 -o mean >  chr12_gc10kb.bed
  
  bedtools map -a sort_hg19_10kb.bed -b  chr13_gc5base.bed  -c 4 -o mean >  chr13_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr14_gc5base.bed  -c 4 -o mean >  chr14_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr15_gc5base.bed  -c 4 -o mean >  chr15_gc10kb.bed
  
  
  bedtools map -a sort_hg19_10kb.bed -b  chr16_gc5base.bed  -c 4 -o mean >  chr16_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr17_gc5base.bed  -c 4 -o mean >  chr17_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr18_gc5base.bed  -c 4 -o mean >  chr18_gc10kb.bed
  
  bedtools map -a sort_hg19_10kb.bed -b  chr19_gc5base.bed  -c 4 -o mean >  chr19_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr20_gc5base.bed  -c 4 -o mean >  chr20_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chr21_gc5base.bed  -c 4 -o mean >  chr21_gc10kb.bed
  
  bedtools map -a sort_hg19_10kb.bed -b  chr22_gc5base.bed  -c 4 -o mean >  chr22_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chrX_gc5base.bed  -c 4 -o mean >  chrX_gc10kb.bed
  bedtools map -a sort_hg19_10kb.bed -b  chrY_gc5base.bed  -c 4 -o mean >  chrY_gc10kb.bed
  
  
  #对 1kb 区间计算GC平均值
  bedtools map -a sort_hg19_1kb.bed -b  chr1_gc5base.bed  -c 4 -o mean >  chr1_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr2_gc5base.bed  -c 4 -o mean >  chr2_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr3_gc5base.bed  -c 4 -o mean >  chr3_gc1kb.bed
  
  bedtools map -a sort_hg19_1kb.bed -b  chr4_gc5base.bed  -c 4 -o mean >  chr4_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr5_gc5base.bed  -c 4 -o mean >  chr5_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr6_gc5base.bed  -c 4 -o mean >  chr6_gc1kb.bed
  
  bedtools map -a sort_hg19_1kb.bed -b  chr7_gc5base.bed  -c 4 -o mean >  chr7_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr8_gc5base.bed  -c 4 -o mean >  chr8_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr9_gc5base.bed  -c 4 -o mean >  chr9_gc1kb.bed
  
  bedtools map -a sort_hg19_1kb.bed -b  chr10_gc5base.bed  -c 4 -o mean >  chr10_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr11_gc5base.bed  -c 4 -o mean >  chr11_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr12_gc5base.bed  -c 4 -o mean >  chr12_gc1kb.bed
  
  bedtools map -a sort_hg19_1kb.bed -b  chr13_gc5base.bed  -c 4 -o mean >  chr13_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr14_gc5base.bed  -c 4 -o mean >  chr14_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr15_gc5base.bed  -c 4 -o mean >  chr15_gc1kb.bed
  
  
  bedtools map -a sort_hg19_1kb.bed -b  chr16_gc5base.bed  -c 4 -o mean >  chr16_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr17_gc5base.bed  -c 4 -o mean >  chr17_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr18_gc5base.bed  -c 4 -o mean >  chr18_gc1kb.bed
  
  
  bedtools map -a sort_hg19_1kb.bed -b  chr19_gc5base.bed  -c 4 -o mean >  chr19_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr20_gc5base.bed  -c 4 -o mean >  chr20_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chr21_gc5base.bed  -c 4 -o mean >  chr21_gc1kb.bed
  
  
  bedtools map -a sort_hg19_1kb.bed -b  chr22_gc5base.bed  -c 4 -o mean >  chr22_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chrX_gc5base.bed  -c 4 -o mean >  chrX_gc1kb.bed
  bedtools map -a sort_hg19_1kb.bed -b  chrY_gc5base.bed  -c 4 -o mean >  chrY_gc1kb.bed
  
  
  #合并每条染色体结果（将上述不同长度区间的文件分别放入不同文件夹）
  files <- list.files()
  read_file <- vector("list", length = 24)
  names(read_file) <- files
  for(i in files)  read_file[[i]] <- fread(i)
  
  for(i in files)
  {
    chr <- strsplit(i, "_gc1kb.bed")[[1]][1]
    df <- read_file[[i]]
    read_file[[i]] <- df[df$V1 == chr,]
  }
  
  names(read_file) <- NULL
  result <- do.call("rbind", read_file)
  
  ###分别输出（）
  
}


##########TFBS##############################
{
  ###其中wgEncodeRegTfbsClusteredWithCellsV3.bed文件为
  ###UCSC下载的全基因组转录因子结合位点数据，论文中有链接地址
  tfbs <- fread("wgEncodeRegTfbsClusteredWithCellsV3.bed")
  bed <- tfbs[,c(1:3,5)]
  write.table(bed,"tfbs.bed", sep = "\t", col.names = F, row.names = F, quote = F)
  
  
  ###human_genome_chromosome_length.tsv文件为从human.hg19.genome文件中获取的1:22+X+Y染色体的长度数据
  sort -k1,1 -k2,2n tfbs.bed > sort_tfbs.bed
  sort -k1,1 -k2,2n human_genome_chromosome_length.tsv > hg19_chr_size.bed
  ###不区分转录因子，计算平均值，同时计算重叠部分的平均值
  bedtools merge -i sort_tfbs.bed -c 4 -o mean > mean_tfbs.bed
  ###将bed格式转换为bigwig格式
  ./bedGraphToBigWig mean_tfbs.bed hg19_chr_size.bed tfbs.bigwig
  
  ###分别计算1Mb,500Kb,...等不同区间长度上的转录因子结合位点数值
  ./bigWigAverageOverBed tfbs.bigwig sort_hg19_mb.bed tfbs_mb.bed
  ./bigWigAverageOverBed tfbs.bigwig sort_hg19_500kb.bed tfbs_500kb.bed
  ./bigWigAverageOverBed tfbs.bigwig sort_hg19_100kb.bed tfbs_100kb.bed
  ./bigWigAverageOverBed tfbs.bigwig sort_hg19_10kb.bed tfbs_10kb.bed
  
}


##########复制时间##########################
{
  bedtools map -a sort_hg19_mb.bed -b sort_average_reptime_14celllines.bed  -c 4 -o mean >  reptime_Mb.bed
  
  bedtools map -a sort_hg19_500kb.bed -b sort_average_reptime_14celllines.bed  -c 4 -o mean >  reptime_500kb.bed
  
  bedtools map -a sort_hg19_100kb.bed -b sort_average_reptime_14celllines.bed  -c 4 -o mean >  reptime_100kb.bed
  
  bedtools map -a sort_hg19_10kb.bed -b sort_average_reptime_14celllines.bed  -c 4 -o mean >  reptime_10kb.bed
}


##########POLII############################
{
  ./bigWigAverageOverBed  wgEncodeSydhTfbsK562Pol2s2StdSig.bigWig  sort_hg19_mb.bed  k562_pol2_Mb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2Etoh01StdSig.bigWig  sort_hg19_mb.bed Mcf10aes_Etoh01_pol2_Mb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2TamStdSig.bigWig  sort_hg19_mb.bed Mcf10aes_Tam_pol2_Mb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsPbdePol2UcdSig.bigWig  sort_hg19_mb.bed Pbde_pol2_Mb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsRajiPol2UcdSig.bigWig  sort_hg19_mb.bed Raji_pol2_Mb.bed
  
  ./bigWigAverageOverBed  wgEncodeSydhTfbsK562Pol2s2StdSig.bigWig  sort_hg19_500kb.bed  k562_pol2_500kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2Etoh01StdSig.bigWig  sort_hg19_500kb.bed Mcf10aes_Etoh01_pol2_500kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2TamStdSig.bigWig  sort_hg19_500kb.bed Mcf10aes_Tam_pol2_500kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsPbdePol2UcdSig.bigWig  sort_hg19_500kb.bed Pbde_pol2_500kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsRajiPol2UcdSig.bigWig  sort_hg19_500kb.bed Raji_pol2_500kb.bed
  
  ./bigWigAverageOverBed  wgEncodeSydhTfbsK562Pol2s2StdSig.bigWig  sort_hg19_100kb.bed  k562_pol2_100kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2Etoh01StdSig.bigWig  sort_hg19_100kb.bed Mcf10aes_Etoh01_pol2_100kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2TamStdSig.bigWig  sort_hg19_100kb.bed Mcf10aes_Tam_pol2_100kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsPbdePol2UcdSig.bigWig  sort_hg19_100kb.bed Pbde_pol2_100kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsRajiPol2UcdSig.bigWig  sort_hg19_100kb.bed Raji_pol2_100kb.bed
  
  ./bigWigAverageOverBed  wgEncodeSydhTfbsK562Pol2s2StdSig.bigWig  sort_hg19_10kb.bed  k562_pol2_10kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2Etoh01StdSig.bigWig  sort_hg19_10kb.bed Mcf10aes_Etoh01_pol2_10kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsMcf10aesPol2TamStdSig.bigWig  sort_hg19_10kb.bed Mcf10aes_Tam_pol2_10kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsPbdePol2UcdSig.bigWig  sort_hg19_10kb.bed Pbde_pol2_10kb.bed
  ./bigWigAverageOverBed  wgEncodeSydhTfbsRajiPol2UcdSig.bigWig  sort_hg19_10kb.bed Raji_pol2_10kb.bed
  
  
  
  setwd("/home/zhangjing/paper/icgc_download/polII/10kb")
  files <- list.files()
  read_file <- vector("list", length = 5)
  names(read_file) <- files
  for(i in files)  read_file[[i]] <- fread(i)
  mean_five <- do.call("cbind", read_file)[,c(5,11,17,23,29)] 
  mean_five <- as.matrix(mean_five)
  mean_value <- apply(mean_five, 1, mean)
  mean_result <- data.frame("num" = read_file[[1]][,1], "length" = read_file[[1]][,2], mean_value)
  write.table(mean_result, "polII_10kb.bed", sep = "\t", row.names = F, quote = F, col.names = F)
  
  
  
  setwd("/home/zhangjing/paper/icgc_download/polII/100kb")
  pol <- fread("polII_100kb.bed")
  mut <- fread("mutation_counts_100kb.tsv")
  result <- merge(mut,pol, by.x = "V4", by.y = "V1")
  result <- result[result$V3.y != 0, ]
  cor(result$V5, result$V3.y)
  
}


##########





