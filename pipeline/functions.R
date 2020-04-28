# 获取区间突变
get_reg_mution = function(mutation_file, region_file, result_file)
{
  temp1 <- paste("sorted_", basename(mutation_file), sep = "")
  temp2 <- paste("sorted_", basename(region_file), sep = "")
  
  cmd1 <- paste("sort -k 1,1 -k 2,2n ",mutation_file," > ",temp1, sep = "")
  cmd2 <- paste("sort -k 1,1 -k 2,2n ",region_file," > ",temp2, sep = "")
  system(cmd1)
  system(cmd2)
  cmd3 <- paste("bedtools intersect -a ",temp1," -b ",temp2," -sorted > ",result_file, sep = "")
  system(cmd3)
  file.remove(temp1, temp2)
}

# 计算区间突变数
get_reg_mution_number = function(mutation_file, region_file, result_file){
  ###将每个突变后面添加数字1
  mut <- fread(mutation_file, data.table = F)
  mut$V4 <- 1
  tempmut = paste0("temp_", basename(mutation_file))
  write.table(mut, tempmut, sep = "\t", row.names = F, col.names = F, quote = F)
  
  ###计算每MB区间内的突变次数
  system(paste0("bedtools map -a ", region_file, " -b ", tempmut, 
         " -c 4 > ", result_file))
  file.remove(tempmut)
  
  ###删除缺失值
  mut_num <- read.table(result_file, stringsAsFactors = F)
  mut_num <- mut_num[mut_num$V5 != ".",]
  mut_num$V5 <- as.numeric(mut_num$V5)
  write.table(mut_num, result_file, sep = "\t", row.names = F, col.names = F, quote = F)
}



#获取CDS的外显子区域
get_cds_region <- function(gtf)
{
  gtf <- fread(gtf, skip = 5)
  
  chr_index <- c(c(1:22),"X","Y")
  gtf <- gtf[gtf$V1 %in% chr_index,]
  gtf_exon <- gtf[gtf$V3 == "exon",]
  
  cds <- gtf_exon[,c(1,4,5)]
  cds$V1 <- paste("chr",as.character(cds$V1),sep = "")
  write.table(cds,"cds_region.bed", sep = "\t", row.names = F, col.names = F, quote = F)
}


###获取bed格式的基因组位点信息
#pos_file format    chr /t start /t end /t ....
#reference_genome_file
pos_fasta <- function(pos_file, ref_genome, result_file)
{
  temp1 <- paste("sorted_", pos_file, sep = "")
  cmd1 <- paste("sort -k 1,1 -k 2,2n ",pos_file," > ",temp1, sep = "")
  system(cmd1)
  cmd2 <- paste("bedtools getfasta -fi ",ref_genome," -bed ",temp1," -tab -bedOut > ",result_file, sep = "")
  system(cmd2)
  file.remove(temp1)
}


###获取bed格式文件的注释文件（bigwig）信息
#pos_file format    chr /t start /t end /t ....
#anno_file format   binary file
pos_anno <- function(pos_file, anno_file, result_file)
{
  temp1 <- paste("sorted_", pos_file, sep = "")
  cmd1 <- paste("sort -k 1,1 -k 2,2n ",pos_file," > ",temp1, sep = "")
  system(cmd1)
  
  bed <- fread(temp1)
  bed <- bed[,1:3]
  bed$index <- 1:dim(bed)[1]
  write.table(bed, "pos_index", sep = "\t", row.names = F, quote = F, col.names = F)
  
  cmd2 <- paste("./bigWigAverageOverBed ",anno_file," pos_index"," index_value", sep = "")
  system(cmd2)
  
  raw_pos <- fread(temp1)
  bigwig_value <- fread("index_value")
  index <- fread("pos_index")
  result <- merge(index, bigwig_value, by.x = "V4", by.y = "V1")
  raw_pos$value <- result$V5
  write.table(raw_pos, result_file, sep = "\t", row.names = F, quote = F, col.names = F)
  
  file.remove(temp1, "pos_index", "index_value")
}


###将hg19基因组分割为不同区间的片段
get_region_bed <- function(region_length)
{
  hg19_chr <- data.table::fread("human.hg19.genome", header = FALSE, data.table = FALSE)
  #hg19_chr <- read.table("human_genome_chromosome_length.tsv", header = F)
  
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
  write.table(result, paste0("hg19_mb_", region_length/1000000L, "M.bed"), sep = "\t", row.names = F, quote = F, col.names = F)
}  


add_coverage = function(gn_region, mut_file, result_file, ref_region="sorted_hg19_1M.bed"){
  tmp_file = paste0("tmp_", basename(gn_region))
  system(paste("bedtools coverage -a", ref_region, "-b", gn_region, ">", tmp_file))
  tmp = fread(tmp_file)
  file.remove(tmp_file)
  mut = fread(mut_file)
  res = merge(mut, tmp, by="V4")
  write.table(res[, c(2:4,1,5,12)], result_file, sep = "\t", row.names = F, col.names = F, quote = F)
}

