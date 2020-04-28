###删除maf函数(将突变位点的1kg AF值大于指定值的位点删除)
#mutation_file format    chr /t start /t end /t ....
rm_maf <- function(threshold, mutation_file)
{
	system("perl extract_1kg_AF.pl")
	system("mv result.tsv 1kg_integrated_phase1_v3_AF.tsv")
	vcf_af <- fread("1kg_integrated_phase1_v3_AF.tsv", sep = "\t", colClasses=list(character= 1), showProgress = F, data.table = F)

	###获取AF大于设定阈值位点
	more_thres <- vcf_af[vcf_af$V6 >= threshold, ]
	pos_index <- paste(paste("chr", more_thres$V1, sep = ""), more_thres$V2, sep = ":")
	
	mut <- fread(mutation_file)
	mut_index <- paste(mut$V1, mut$V3, sep = ":")
	common_pos <- intersect(pos_index, mut_index)
	rm_snp_mut <- mut[!mut_index %in% common_pos,]
	write.table(rm_snp_mut, paste("rmSNP", mutation_file, sep = "_"), sep = "\t", row.names = F, col.names = F,quote = F)
	
	file.remove("1kg_integrated_phase1_v3_AF.tsv")
}


###求两个区间的交集
# file1, file2 format   chr /t start /t end /t ....  
get_intersect_region <- function(file1, file2, file_result)
{
	temp1 <- paste("sorted_", file1, sep = "")
	temp2 <- paste("sorted_", file2, sep = "")
	cmd1 <- paste("sort -k 1,1 -k 2,2n ",file1," > ",temp1, sep = "")
	cmd2 <- paste("sort -k 1,1 -k 2,2n ",file2," > ",temp2, sep = "")
	system(cmd1)
	system(cmd2)
	
	cmd3 <- paste("bedtools intersect -a ",temp1," -b ",temp2," -sorted > ",file_result, sep = "")
	system(cmd3)
	file.remove(temp1, temp2)	
}



###获取指定区间的突变
#mutation_file format    chr /t start /t end /t ....
#region_file format    chr /t start /t end /t ....
region_mut <- function(mutation_file, region_file, result_file)
{
	temp1 <- paste("sorted_", mutation_file, sep = "")
	temp2 <- paste("sorted_", region_file, sep = "")
	
	cmd1 <- paste("sort -k 1,1 -k 2,2n ",mutation_file," > ",temp1, sep = "")
	cmd2 <- paste("sort -k 1,1 -k 2,2n ",region_file," > ",temp2, sep = "")
	system(cmd1)
	system(cmd2)
	cmd3 <- paste("bedtools intersect -a ",temp1," -b ",temp2," -sorted > ",result_file, sep = "")
	system(cmd3)
	file.remove(temp1, temp2)
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




###获取bed格式文件的注释文件（bed）信息
#pos_file format    chr /t start /t end /t ....
#anno_file format   chr /t start /t end /t ....
pos_anno <- function(pos_file, anno_file, result_file)
{
    temp1 <- paste("sorted_", pos_file, sep = "")
	temp2 <- paste("sorted_", anno_file, sep = "")
	cmd1 <- paste("sort -k 1,1 -k 2,2n ",pos_file," > ",temp1, sep = "")
	cmd2 <- paste("sort -k 1,1 -k 2,2n ",anno_file," > ",temp2, sep = "")
	system(cmd1)
	system(cmd2)
	cmd3 <- paste("bedtools map -a ",temp1," -b ",temp2," -c 4 -o mean > ",result_file, sep = "")
	system(cmd3)
	file.remove(temp1, temp2)
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










###筛选需要的列
raw_mut_file <- "simple_somatic_mutation.open.tsv.gz"
vec_string1 <- c("zcat ", raw_mut_file, " | ", 
				 "awk -F'\t' -vOFS='\t' '{print $1,$2,$3,$9,$10,$11,$12,$14,$15,$16,$17,$22,$26,$29,$30}' ",
				 "> extract_simple_somatic_mutation.open.tsv")
cmd1 <- paste(vec_string1, collapse = "")
system(cmd1)






###删除超过500,000突变的病例
threshold <- 500000

#Remove ultra-mutated donors' mutation data
system("awk -F'\t' -vOFS='\t' '{print $1,$2}' extract_simple_somatic_mutation.open.tsv > mutId_donor.tsv")
id_donor <- fread("mutId_donor.tsv", showProgress = F, data.table = F)
id_donor <- unique(id_donor)
donor_mutNum <- as.data.frame(table(id_donor$icgc_donor_id))
#get ultra-mutated donors' ID
ultra_mut_donor <- donor_mutNum[donor_mutNum[,2] > threshold, ]
ultra_donor <- as.character(ultra_mut_donor[,1])
	
vec_string2 <- c("sed '",
				 paste(paste("/",ultra_donor, "/d", sep = ""), collapse = ";"),					 
				 "' extract_simple_somatic_mutation.open.tsv > rmHyperMut_ext_ssm.open.tsv")
cmd2 <- paste(vec_string2, collapse = "")				
system(cmd2)
	

	
###删除重复的行
system("awk -F'\t' -vOFS='\t' '{print $1,$2,$4,$5,$6,$7,$8,$9,$10,$11}' rmHyperMut_ext_ssm.open.tsv > sub_rmHyperMut_ext_ssm.tsv")  

rm_redundant <- fread("sub_rmHyperMut_ext_ssm.tsv", sep = "\t", data.table = F)
all_point_mut <- unique(rm_redundant)
names(all_point_mut) <- c("mut_id","donor","chr","start","end","ref", "mut")
all_point_mut$chr <- paste("chr",all_point_mut$chr, sep = "")
fwrite(all_point_mut,"预处理所有突变.tsv", sep = "\t", row.names = F, quote = F)



###选取 "预处理所有突变.tsv" 文件中的单位点突变
mut <- fread("预处理所有突变.tsv", data.table = F)
pos_mut <- mut[mut$mut == "single base substitution",]



###以bed格式输出所有单碱基突变
mut_bed <- pos_mut[,c(3:5,2)]
mut_bed$start <- mut_bed$start - 1
write.table(mut_bed, "single_mut.bed", sep = "\t", row.names = F, quote = F, col.names = F)



###删除所有突变中的snp大于0.05的突变
rm_maf(0.05, "single_mut.bed")
file.remove("single_mut.bed")




#获取高拼接性区域
get_low_map_region <- function()
{
	thres <- 1
	system("./bigWigToBedGraph wgEncodeCrgMapabilityAlign24mer.bigWig wgEncodeCrgMapabilityAlign24mer.bedGraph")
	mer100 <- fread("wgEncodeCrgMapabilityAlign24mer.bedGraph")
	mer100_thres <- mer100[mer100$V4 == thres,]
	write.table(mer100_thres,"high_map.bed", sep = "\t", row.names = F, quote = F, col.names = F)
	
	file.remove("wgEncodeCrgMapabilityAlign24mer.bedGraph")
}






#获取CDS的外显子区域
get_cds_region <- function()
{
	gtf <- fread("Homo_sapiens.GRCh37.75.gtf", skip = 5)
	
	chr_index <- c(c(1:22),"X","Y")
	gtf <- gtf[gtf$V1 %in% chr_index,]
	gtf_exon <- gtf[gtf$V3 == "exon",]
	
	cds <- gtf_exon[,c(1,4,5)]
	cds$V1 <- paste("chr",as.character(cds$V1),sep = "")
	write.table(cds,"cds_region.bed", sep = "\t", row.names = F, col.names = F, quote = F)
}




get_hyperimmu_region <- function()
{
	###从Ensmbles gtf 文件 Homo_sapiens.GRCh37.75.gtf 中获取免疫高变区区间
	immu_id <- c('IG_C_gene',
				'IG_D_gene',
				'IG_J_gene',
				'IG_LV_gene',
				'IG_V_gene',
				'TR_C_gene',
				'TR_J_gene',
				'TR_V_gene',
				'TR_D_gene',
				'IG_pseudogene',
				'IG_C_pseudogene',
				'IG_J_pseudogene',
				'IG_V_pseudogene',
				'TR_V_pseudogene',
				'TR_J_pseudogene')
				
	gtf <- fread("Homo_sapiens.GRCh37.75.gtf", skip = 5)

	
	chr_index <- c(c(1:22),"X","Y")
	gtf <- gtf[gtf$V1 %in% chr_index,]
	
	df_immu_region <- gtf[gtf$V2 %in% immu_id,]
	bed_immu_region <- data.frame(df_immu_region$V1,df_immu_region$V4, df_immu_region$V5)
	names(bed_immu_region) <- c("chr","start","end")
	bed_immu_region$chr <- paste("chr",as.character(bed_immu_region$chr),sep = "")
	write.table(bed_immu_region,"immu_region_from_gtf.bed", sep = "\t", row.names = F, col.names = F,quote = F)
	
	###将免疫高变区merge
	system("sort -k 1,1 -k 2,2n immu_region_from_gtf.bed > sort_chr_immu_region_from_gtf.bed")
	system("bedtools merge -i sort_chr_immu_region_from_gtf.bed > merged_sort_chr_immu_region_from_gtf.bed")

	
	###将免疫高变区前后扩增ext KB
	ext = 50000
	merged <- fread("merged_sort_chr_immu_region_from_gtf.bed")
	merged$V2 <- merged$V2 - ext
	merged$V3 <- merged$V3 + ext
	write.table(merged,"extended_merged_sort_chr_immu_region_from_gtf.bed", sep = "\t", row.names = F ,col.names = F, quote = F)

	###再次merge
	system("sort -k 1,1 -k 2,2n extended_merged_sort_chr_immu_region_from_gtf.bed > sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
	system("bedtools merge -i sorted_extended_merged_sort_chr_immu_region_from_gtf.bed > merged_sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
	
	file.remove("sort_chr_immu_region_from_gtf.bed",
				"merged_sort_chr_immu_region_from_gtf.bed",
				"immu_region_from_gtf.bed",
				"sorted_extended_merged_sort_chr_immu_region_from_gtf.bed",
				"extended_merged_sort_chr_immu_region_from_gtf.bed")
}



###获取上面两个区域的并集区间
get_mask_region <- function()
{
	cds <- fread("cds_region.bed")
   #low_map <- fread("low_map.bed")[,1:3]
	immu <- fread("merged_sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
	result <- rbind(cds,immu)
	write.table(result, "rbind_three.bed", sep = "\t", row.names = F, col.names = F, quote = F)	
	
	system("sort -k 1,1 -k 2,2n rbind_three.bed > sorted_rbind_three.bed")
	system("bedtools merge -i sorted_rbind_three.bed > merged_sorted_rbind_three.bed")
	
	file.remove("rbind_three.bed","sorted_rbind_three.bed")
}





###获取基因组上除去上述三个区间的剩下的非编码区间
get_nonmask_region <- function()
{
	system("sort -k 1,1 -k 2,2n merged_sorted_rbind_three.bed > sorted_merged_sorted_rbind_three.bed")
	system("sort -k 1,1 -k 2,2n human.hg19.genome > sorted_human.hg19.genome")
	######获取全基因组mask区间的补集######
	system("bedtools complement -i sorted_merged_sorted_rbind_three.bed -g sorted_human.hg19.genome > hg19_nonmask_region.bed")
	######获取补集区间与高map区间的交集
	get_intersect_region("hg19_nonmask_region.bed","high_map.bed","hg19_map_nonmask_region.bed")
	
	
	file.remove("sorted_merged_sorted_rbind_three.bed","sorted_human.hg19.genome", "hg19_nonmask_region.bed")
}



###获取非编码区的所有体细胞突变
options(scipen = 30)
write.table(fread("rmSNP_single_mut.bed"),"rmSNP_single_mut.bed", sep = "\t", row.names = F, col.names = F, quote = F)
region_mut("rmSNP_single_mut.bed","hg19_map_nonmask_region.bed","noncoding_mut")


###随机获取基因组上50000000个位点
system("bedtools random -g human.hg19.genome -l 1 -n 50000000 > hg19_random_site.bed")
###获取非编码区的所有随机突变
region_mut("hg19_random_site.bed","hg19_map_nonmask_region.bed","noncoding_random_site.bed")
file.remove("hg19_random_site.bed")





###从非编码随机位点中将体细胞突变位点删除并将突细胞突变位点和随机位点合并
noncod_mut <- fread("noncoding_mut")
noncod_mut <- noncod_mut[noncod_mut$V1 != "chrY",]
noncod_site <- fread("noncoding_random_site.bed")[,1:3]
noncod_site$V4 <- "mock"

mut_index <- paste(noncod_mut$V1, noncod_mut$V3, sep = ":")
site_index <- paste(noncod_site$V1, noncod_site$V3, sep = ":")
common_pos <- intersect(mut_index, site_index)
noncod_site_rm_mut <- noncod_site[!site_index %in% common_pos, ]
noncod_site_rm_mut$index <- 0
noncod_mut$index <- 1
all_site <- rbind(noncod_mut, noncod_site_rm_mut)
write.table(all_site, "noncoding_mut_site.bed", sep = "\t", row.names = F, quote = F, col.names = F)








###获取 noncoding_mut_site.bed 文件的三碱基
single_site <- fread("noncoding_mut_site.bed")
single_site$V2 <- single_site$V2 - 1
single_site$V3 <- single_site$V3 + 1
write.table(single_site, "three_base.bed", sep = "\t", row.names = F, col.names = F, quote = F)

pos_fasta("three_base.bed","hg19.fasta","noncoding_three_base.bed")
file.remove("three_base.bed")


three_site <- fread("noncoding_three_base.bed")
three_site$V2 <- three_site$V2 + 1
three_site$V3 <- three_site$V3 - 1
three_site$V6 <- toupper(three_site$V6)

base_type <- fread("base.tsv")
merge_result <- merge(three_site, base_type, by.x = "V6", by.y = "V1")
result <- merge_result[,c(2:6,1,7)]
write.table(result, "noncoding_mut_site_3basetype.bed", sep = "\t", row.names = F, col.names = F, quote = F)






###获取复制时间值(chr start end value)
get_reptime <- function(path_to_rep_file)
{
	setwd(path_to_rep_file)
	system("./bigWigToWig wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig BJ_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig GM06990_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig GM12801_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig GM12812_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig GM12813_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig GM12878_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig HelaS3_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig Hepg2_reptime_wavesignal.wig")		
	system("./bigWigToWig wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig Huvec_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig Imr90_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig K562_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig Mcf7_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig Nhek_reptime_wavesignal.wig")
	system("./bigWigToWig wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig Sknsh_reptime_wavesignal.wig")

	dir.create("wigfile")
	system("mv *.wig ./wigfile")
	dir <- paste(path_to_rep_file, "/wigfile", sep = "")
	setwd(dir)
		
	file_names <- list.files()
	file_anno <- vector("list", length = length(file_names))
	names(file_anno) <- file_names
	
	#function to get comment lines for each wig files
	get_anno_line <- function(file_name)
	{
		con_file <- file(file_name, "r")
		read_line <- readLines(con_file)
		anno_line <- read_line[grep("fixedStep", read_line)]
		close(con_file)
		anno_line
	}
	
	#check if the divided reptime regions of all cell lines are identical 
	for(i in file_names)  file_anno[[i]] <- get_anno_line(i)
	is_identical <- T
	for(i in 2:length(file_names)) 
	{
		if(!identical(file_anno[[i]] , file_anno[[1]]))  
		{
			is_identical <- F
			break
		}
	}
	if(is_identical == T){
		print("ALL replication time files are the same region")
	}else{
		print("ALL replication time files are not he same region")
	}
	
	
	#function to get the reptime for each wig file ,the output is bed format
	get_bed_format <- function(wig_file)
	{
		read_file <- file(wig_file, "r")
		read_lines <- readLines(read_file)
		close(read_file)
		
		#get the index of row  of all "fixedStep" lines
		fixedline_index <- grep("fixedStep", read_lines)
		
		interval_end <- fixedline_index[2:length(fixedline_index)] - 1
		interval_end <- c(interval_end, length(read_lines))
		
		fixedline <- read_lines[fixedline_index]
		bed_df_interval <- vector("list", length = length(fixedline_index))
		names(bed_df_interval) <- fixedline

		for(i in 1:length(fixedline))
		{
			chr = strsplit(strsplit(fixedline[i],"chrom=")[[1]][2]," ")[[1]][1]
			start_pos <- strsplit(strsplit(fixedline[i],"start=")[[1]][2]," ")[[1]][1]
			interval_val <- read_lines[(fixedline_index[i]+1):interval_end[i]]
			row_num <- length(interval_val)
			bed_start_pos <- seq(as.numeric(start_pos),by = 1000, length.out = row_num)
			bed_end_pos <- bed_start_pos + 1000
			bed_df_interval[[fixedline[i]]] <- data.frame(chr, bed_start_pos, bed_end_pos, interval_val)
		}
		
		names(bed_df_interval) <- NULL
		unlist_bed <- list(rbindlist(bed_df_interval))
		names(unlist_bed) <- wig_file
		unlist_bed
	}
	
	#get bed format reptime for all wig files
	bed_fomart_list <- lapply(file_names, get_bed_format)
	
	#Check if all elments in bed format list are the same region
	all_identical <- T
	for(i in 2:length(bed_fomart_list))
	{
		if(!identical(bed_fomart_list[[i-1]][[1]][,1:3] , bed_fomart_list[[i]][[1]][,1:3])) 
		{
			all_identical <- F
			break
		}
	}
	if(all_identical == T){
		print("ALL file have the same format and region")
	}else{
		print("ALL file have not the same format and region")
	}

	#calculate average reptime for all celllines
	val_mat <- matrix(nrow = dim(bed_fomart_list[[1]][[1]])[1], ncol = length(bed_fomart_list))
	for(i in 1:length(bed_fomart_list))  
		val_mat[,i] <- as.numeric(as.character((bed_fomart_list[[i]][[1]]$interval_val)))
	mean_vec <- apply(val_mat,1, mean)

	###transform to bed format
	average_rep_result <- data.frame( bed_fomart_list[[1]][[1]][,1:3], mean_vec)
	names(average_rep_result) <- c("chr","start", "end", "rep_val")
	write.table(average_rep_result, "average_reptime_14celllines.bed", sep = "\t", row.names = F, quote = F)
}




###获取 noncoding_mut_site_3basetype.bed 文件的复制时间
pos_anno("noncoding_mut_site_3basetype.bed","average_reptime_14celllines.bed","noncoding_mut_site_3basetype_reptime.bed")
	
rep_value <- fread("noncoding_mut_site_3basetype_reptime.bed")
rep_value <- rep_value[rep_value$V8 != ".",]
write.table(rep_value, "noncoding_mut_site_3basetype_reptime.bed", sep = "\t", row.names = F, quote = F, col.names = F)


###获取 noncoding_mut_site_3basetype_reptime.bed 文件的tfbs
tfbs <- fread("wgEncodeRegTfbsClusteredWithCellsV3.bed")
bed <- tfbs[,c(1:3,5)]
write.table(bed,"TFBS_value.bed", sep = "\t", col.names = F, row.names = F, quote = F)

pos_anno("noncoding_mut_site_3basetype_reptime.bed","TFBS_value.bed","noncoding_mut_site_3basetype_reptime_tfbs.bed")
tfbs_value <- fread("noncoding_mut_site_3basetype_reptime_tfbs.bed")
tfbs_value$V9[tfbs_value$V9 == "."] <- 0
write.table(tfbs_value, "noncoding_mut_site_3basetype_reptime_tfbs.bed", sep = "\t", row.names = F, quote = F, col.names = F)



###获取 noncoding_mut_site_3basetype_reptime_tfbs.bed 文件位点的保守性值
pos_anno("noncoding_mut_site_3basetype_reptime_tfbs.bed","hg19.100way.phastCons.bw","noncoding_mut_site_3basetype_reptime_tfbs_cons.bed")


###获取 noncoding_mut_site_3basetype_reptime_tfbs.bed 文件位点的GC含量（1kb）
system("sort -k 1,1 -k 2,2n hg19_gc1kb.bed > sorted_hg19_gc1kb.bed")
gc <- fread("sorted_hg19_gc1kb.bed")
gc <- gc[gc$V5!= ".",]
write.table(gc, "sorted_hg19_gc1kb.bed", sep = "\t", row.names = F, col.names = F, quote = F)
system("bedtools map -a noncoding_mut_site_3basetype_reptime_tfbs_cons.bed -b sorted_hg19_gc1kb.bed -c 5 -o mean > noncoding_mut_site_3basetype_reptime_tfbs_cons_gc.bed")

file.remove("sorted_hg19_gc1kb.bed")



###获取启动子区域并判断是否在启动子区
	###选取-2500 到 +500 作为启动子区域
	up_site = 2500
	down_site = 500
	##过滤数据
	gtf_ensmble <- fread("Homo_sapiens.GRCh37.75.gtf", sep = "\t", colClasses=list(character= 1))
	names(gtf_ensmble)[1] <- "V1"
	ensmble_pro_cod <- gtf_ensmble[gtf_ensmble$V2 == "protein_coding",]            #############选择仅仅是蛋白质编码的转录本
	transcript_pro_cod_ensmble <- ensmble_pro_cod[ensmble_pro_cod$V3 == "transcript",]          ##########从转录本中选取总的转录本区间
	options(scipen=30)
	chr_index <- c(c(1:22),"X","Y")
	chr_trans_pro <- transcript_pro_cod_ensmble[transcript_pro_cod_ensmble$V1 %in% chr_index, ]       #########仅仅挑选chr1:22 + X + Y 染色体
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

	system("sort -k 1,1 -k 2,2n promoter_region.bed > sorted_promoter_region.bed")
	system("bedtools map -a noncoding_mut_site_3basetype_reptime_tfbs_cons_gc.bed -b sorted_promoter_region.bed -c 4 -o mean > noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter.bed")
	
	pro <- fread("noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter.bed")
	pro$V12[pro$V12 == "."] <- "N"
	pro$V12[pro$V12 == 1] <- "Y"
	write.table(pro, "noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter.bed", sep = "\t", row.names = F, col.names = F, quote = F)

	
	file.remove("promoter_region.bed","sorted_promoter_region.bed")





	
###判断突变是不是在CpG岛
	cpg <- fread("CpG_island_UCSC.tsv")
	chr <- c(1:22,"X", "Y")
	chr <- paste("chr", chr, sep= "")
	cpg <- cpg[cpg$chrom %in% chr,]
	bed <- cpg[,2:4]
	bed$index <- 1
	write.table(bed,"cpg.bed", sep = "\t", row.names = F, quote = F, col.names = F)
	
	system("sort -k 1,1 -k 2,2n cpg.bed > sorted_cpg.bed")
	system("bedtools map -a noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter.bed -b sorted_cpg.bed -c 4 -o mean > noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter_cpg.bed")
	
	cpg <- fread("noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter_cpg.bed")
	cpg$V13[cpg$V13 == "."] <- "N"
	cpg$V13[cpg$V13 == 1] <- "Y"
	write.table(cpg, "noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter_cpg.bed", sep = "\t", row.names = F, quote = F, col.names = F)

	file.remove("cpg.bed","sorted_cpg.bed")
	
	

	
system("mv noncoding_mut_site_3basetype_reptime_tfbs_cons_gc_promoter_cpg.bed  genetic.bed")

	
	
	
	
###获取不同的组织的表观遗传注释
broadpeak2bed <- function()
{
  files <- list.files()
  files <- files[grep(".broadPeak", files)]
  
  read_file <- vector("list", length = length(files))
  names(read_file) <- files
  for(i in files)  read_file[[i]] <- fread(i)
  for(i in files)  
  {
    bed <- read_file[[i]][,c(1:3,7)]
    write.table(bed, paste(i,"bed", sep = "."), sep = "\t", row.names = F, col.names = F, quote = F)
  }
}

pos_anno <- function(pos_file, anno_file, result_file)
{
	temp2 <- paste("sorted_", anno_file, sep = "")
	cmd2 <- paste("sort -k 1,1 -k 2,2n ",anno_file," > ",temp2, sep = "")
	system(cmd2)
	cmd3 <- paste("bedtools map -a ",pos_file," -b ",temp2," -c 4 -o mean > ",result_file, sep = "")
	system(cmd3)
	file.remove(temp2)
}








library(stringr)

setwd("/home/zhangjing/roadmap/blood")
setwd("/home/zhangjing/roadmap/breast")
setwd("/home/zhangjing/roadmap/Esophagus")
setwd("/home/zhangjing/roadmap/Kidney")
setwd("/home/zhangjing/roadmap/liver")
setwd("/home/zhangjing/roadmap/lung")
setwd("/home/zhangjing/roadmap/Ovary")
setwd("/home/zhangjing/roadmap/Pancreas")



broadpeak2bed()
system("mkdir bed")
system("mv *.broadPeak.bed bed")
setwd("bed")


file_names <- list.files()
anno <- str_split_fixed(file_names, "[-|.]",4)[,2]
anno_file <- vector(length = length(file_names) + 1)
anno_file[1] <- "genetic.bed"
for(i in 1:length(file_names))  anno_file[i + 1] <- paste(anno_file[i], anno[i], sep = ".")


****mv genetic.bed to current direct


for(i in 1:length(file_names))    
{
	pos_anno(anno_file[i], file_names[i], anno_file[i+1])
}



result_file <- fread(anno_file[length(file_names) + 1])
names(result_file)[14: (13 + length(file_names))] <- anno
result_file[result_file == "."]  <- 0
write.table(result_file,"final_result.tsv", sep = "\t", row.names = F, quote = F)
system("mv genetic.bed ..")
file.remove(anno_file)







######################################################## 多个数据

###melanoma
setwd("/home/zhangjing/roadmap/E059_Melanoma")


broadpeak2bed()
system("mkdir bed")
system("mv *.broadPeak.bed bed")
setwd("bed")



file_names <- list.files()
anno <- str_split_fixed(file_names, ".broadPeak.bed",4)[,1]
anno_file <- vector(length = length(file_names) + 1)
anno_file[1] <- "genetic.bed"
for(i in 1:length(file_names))  anno_file[i + 1] <- paste(anno_file[i], anno[i], sep = ".")


for(i in 1:length(file_names))    
{
	pos_anno(anno_file[i], file_names[i], anno_file[i+1])
}



result_file <- fread(anno_file[length(file_names) + 1])
names(result_file)[14: (13 + length(file_names))] <- anno
result_file[result_file == "."]  <- 0
write.table(result_file,"final_result.tsv", sep = "\t", row.names = F, quote = F)
system("mv genetic.bed ..")
file.remove(anno_file)




setwd("/home/zhangjing/roadmap/E061_Melanoma")


broadpeak2bed()
system("mkdir bed")
system("mv *.broadPeak.bed bed")
setwd("bed")



file_names <- list.files()
anno <- str_split_fixed(file_names, ".broadPeak.bed",4)[,1]
anno_file <- vector(length = length(file_names) + 1)
anno_file[1] <- "genetic.bed"
for(i in 1:length(file_names))  anno_file[i + 1] <- paste(anno_file[i], anno[i], sep = ".")


for(i in 1:length(file_names))    
{
	pos_anno(anno_file[i], file_names[i], anno_file[i+1])
}



result_file <- fread(anno_file[length(file_names) + 1])
names(result_file)[14: (13 + length(file_names))] <- anno
result_file[result_file == "."]  <- 0
write.table(result_file,"final_result.tsv", sep = "\t", row.names = F, quote = F)
system("mv genetic.bed ..")
file.remove(anno_file)








install.packages("pROC")
install.packages("caret")


#十折交叉验证和ROC曲线
###肺癌
{
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/lung/bed")
system("mv final_result.tsv cal_model_lung_final_result.tsv")



project <- c("LUSC-CN","LUSC-KR","mock")
lung_mut <- logit_form(project, "cal_model_lung_final_result.tsv")


###多重共线性检验
mt <- lung_mut[,c(8:11,14:19)] 
x <- cor(mt)         
kappa(x, exact = TRUE)         # 8.947881                   
mt <- lung_mut[,c(8:11)]
x <- cor(mt)
kappa(x, exact = TRUE)         # 3.16644        
mt <- lung_mut[,c(14:19)]
x <- cor(mt)
kappa(x, exact = TRUE)		   # 5.184799


#检验数值类型
is.character(lung_mut[,7])
is.numeric(as.matrix(lung_mut[,c(8:11,14:19)]))


#整理格式及结果计算
lung_form <- lung_mut[,c(5,7:19)]
cross_val(lung_form,"lung_logit.RData")                       # 0.6850658           
}                                  



##食道癌
{
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/Esophagus/bed")
system("mv final_result.tsv cal_model_esophagus_final_result.tsv")

project <- c("ESAD-UK","ESCA-CN","mock")
esophagus_mut <- logit_form(project, "cal_model_esophagus_final_result.tsv")


###多重共线性检验
mt <- esophagus_mut[,c(8:11,14:19)] 
x <- cor(mt)
kappa(x, exact = TRUE)            # 8.529881           
mt <- esophagus_mut[,c(8:11)]
x <- cor(mt)
kappa(x, exact = TRUE)            # 3.178326      
mt <- esophagus_mut[,c(14:19)]
x <- cor(mt)
kappa(x, exact = TRUE)			  # 4.565268



#检验数值类型
is.character(esophagus_mut[,7])
is.numeric(as.matrix(esophagus_mut[,c(8:11,14:19)]))

#整理格式
esophagus_form <- esophagus_mut[,c(5,7:19)]
cross_val(esophagus_form,"esophagus_logit.RData")
}



####肝癌
{
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/liver/bed")
system("mv final_result.tsv cal_model_liver_final_result.tsv")


project <- c("LIAD-FR", "LICA-CN", "LICA-FR", "LINC-JP", "LIRI-JP","mock")
liver_mut <- logit_form(project, "cal_model_liver_final_result.tsv")


###多重共线性检验
mult_validate(liver_mut)	 # 8.991961		

#检验数据类型
type_validate(liver_mut)   

#整理格式
liver_form <-format_trans(liver_mut)
cross_val(liver_form,"liver_logit.RData")   
}



{#######乳腺癌
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/breast/bed")
system("mv final_result.tsv cal_model_breast_final_result.tsv")

project <- c("BRCA-EU","BRCA-FR","BRCA-US","mock")
breast_mut <- logit_form(project, "cal_model_breast_final_result.tsv")

#添加dnase数据
a <- fread("E028-DNase.hotspot.broad.bed")[,-4]
write.table(a,"E028-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)


add_res <- add_dnase(breast_mut, "E028-DNase.hotspot.broad.bed")
breast_mut <- add_res

###多重共线性检验
mult_validate(breast_mut)

#检验数据类型
type_validate(breast_mut)


#整理格式
breast_form <-format_trans(breast_mut)
cross_val(breast_form,"breast_logit.RData")
}



## 胰腺癌
{
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/Pancreas/bed")
system("mv final_result.tsv cal_model_pancreas_final_result.tsv")
	

project <- c("PACA-AU","PACA-CA","PAEN-AU","PAEN-IT","mock")
pancreas_mut <- logit_form(project, "cal_model_pancreas_final_result.tsv")

#添加dnase数据
a <- fread("E098-DNase.hotspot.broad.bed")[,-4]
write.table(a,"E098-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)

add_res <- add_dnase(pancreas_mut, "E098-DNase.hotspot.broad.bed")
pancreas_mut <- add_res

###多重共线性检验
mult_validate(pancreas_mut)
				
#检验数据类型
type_validate(pancreas_mut)

#整理格式
pancreas_form <-format_trans(pancreas_mut)
cross_val(pancreas_form,"pancreas_logit.RData")
}



###肾癌
{
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/Kidney/bed")
system("mv final_result.tsv cal_model_kidney_final_result.tsv")


project <- c("RECA-EU","mock")
kidney_mut <- logit_form(project, "cal_model_kidney_final_result.tsv")

#添加dnase数据
a <- fread("E086-DNase.hotspot.broad.bed")[,-4]
write.table(a,"E086-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)

add_res <- add_dnase(kidney_mut, "E086-DNase.hotspot.broad.bed")
kidney_mut <- add_res

###多重共线性检验
mult_validate(kidney_mut)			

#检验数据类型
type_validate(kidney_mut)

#整理格式
kidney_form <- format_trans(kidney_mut)
# cross_val(kidney_form,"kidney_logit.RData")    #AUC 0.626 不合格,优化

##调整    res ~ nbase + reptime + conser + cpg + H3K27me3 + H3K36me3 + H3K9me3  
kidney_form1 <- kidney_form[,c(1,2,3,5,8,9,10,14)]     ###使用该结果
cross_val(kidney_form1,"kidney1_logit.RData") 
}



{####血癌
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/blood/bed")
system("mv final_result.tsv cal_model_blood_final_result.tsv")


project <- c("ALL-US","CLLE-ES","MALY-DE","NKTL-SG","mock")
blood_mut <- logit_form(project, "cal_model_blood_final_result.tsv")


###多重共线性检验
mult_validate(blood_mut)			

#检验数据类型
type_validate(blood_mut)

#整理格式
blood_form <-format_trans(blood_mut)
cross_val(blood_form,"blood_logit.RData")  
}



###肝癌
{
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/liver/bed")
system("mv final_result.tsv cal_model_liver_final_result.tsv")


project <- c("LIAD-FR", "LICA-CN", "LICA-FR", "LINC-JP", "LIRI-JP","mock")
liver_mut <- logit_form(project, "cal_model_liver_final_result.tsv")


###多重共线性检验
mult_validate(liver_mut)			

#检验数据类型
type_validate(liver_mut)

#整理格式
liver_form <-format_trans(liver_mut)
cross_val(liver_form,"liver_logit.RData")  

}



###卵巢癌
{
source("/home/zhangjing/paper/lasso/logit_function20181112.R")
setwd("/home/zhangjing/roadmap/Ovary/bed")
system("mv final_result.tsv cal_model_ovary_final_result.tsv")

project <- c("OV-AU","mock")
ovary_mut <- logit_form(project, "cal_model_ovary_final_result.tsv")


a <- fread("E097-DNase.hotspot.broad.bed")[,-4]
write.table(a,"E097-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)

add_res <- add_dnase(ovary_mut, "E097-DNase.hotspot.broad.bed")
ovary_mut <- add_res


###多重共线性检验
mult_validate(ovary_mut)			

#检验数据类型
type_validate(ovary_mut)

#整理格式
ovary_form <-format_trans(ovary_mut)
#cross_val(ovary_form,"ovary_logit.RData")   
###更改
ovary_form1 <- ovary_form[,-c(6:8,15)]
cross_val(ovary_form1,"ovary1_logit.RData")  







}
 

 
####黑色素瘤
{
	###获取E059和E061的平均值
	setwd("/home/zhangjing/roadmap/E059_Melanoma/bed")
	system("mv final_result.tsv E059_final_result.tsv")
	setwd("/home/zhangjing/roadmap/E061_Melanoma/bed")
	system("mv final_result.tsv E061_final_result.tsv")
	
	df59 <- fread("/home/zhangjing/roadmap/E059_Melanoma/bed/E059_final_result.tsv", data.table = F)
	df61 <- fread("/home/zhangjing/roadmap/E061_Melanoma/bed/E061_final_result.tsv", data.table = F)
	e059 <- as.matrix(df59[,14:19])
	e061 <- as.matrix(df61[,14:19])
	
	mean_value <- (e059 + e061)/2
	mean_value <- as.data.frame(mean_value)
	
	names(mean_value) <- c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3")
	res <- data.frame(df59[,1:13], mean_value)
	write.table(res, "final_result.tsv", sep = "\t", row.names = F, quote = F)

	
	
	source("/home/zhangjing/paper/lasso/logit_function20181112.R")
	setwd("/home/zhangjing/roadmap/E061_Melanoma/bed")
	system("mv final_result.tsv cal_model_mela_final_result.tsv")
	
	project <- c("MELA-AU","SKCA-BR","SKCM-US","mock")
	mela_mut <- logit_form(project, "cal_model_mela_final_result.tsv")

	#添加dnase数据
	a <- fread("E059-DNase.hotspot.broad.bed")[,-4]
	write.table(a,"E059-DNase.hotspot.broad.bed", sep = "\t", row.names = F, quote = F, col.names = F)

	add_res <- add_dnase(mela_mut, "E059-DNase.hotspot.broad.bed")
	mela_mut <- add_res

	###多重共线性检验
	mult_validate(mela_mut)			

	#检验数据类型
	type_validate(mela_mut)

	#整理格式
	mela_form <-format_trans(mela_mut)
	cross_val(mela_form,"mela_logit.RData")   	
}






###取所有的突变并删除不符合条件的位点
setwd("/home/zhangjing/predict_prob")
{###选取 "预处理所有突变.tsv" 文件中的单位点突变
mut <- fread("预处理所有突变.tsv", data.table = F)
pos_mut <- mut[mut$mut == "single base substitution",]


###以bed格式输出所有单碱基突变
mut_bed <- pos_mut[,c(3:5,2)]
mut_bed$start <- mut_bed$start - 1
write.table(mut_bed, "single_mut.bed", sep = "\t", row.names = F, quote = F, col.names = F)
}


###删除maf0.05
rm_maf(0.05, "single_mut.bed")



###获取免疫高变区，并扩增指定bp
{
	###从Ensmbles gtf 文件 Homo_sapiens.GRCh37.75.gtf 中获取免疫高变区区间
	immu_id <- c('IG_C_gene',
				'IG_D_gene',
				'IG_J_gene',
				'IG_LV_gene',
				'IG_V_gene',
				'TR_C_gene',
				'TR_J_gene',
				'TR_V_gene',
				'TR_D_gene',
				'IG_pseudogene',
				'IG_C_pseudogene',
				'IG_J_pseudogene',
				'IG_V_pseudogene',
				'TR_V_pseudogene',
				'TR_J_pseudogene')
				
	gtf <- fread("Homo_sapiens.GRCh37.75.gtf", sep = "\t", colClasses=list(character= 1))
	names(gtf)[1] <- "V1"
	
	chr_index <- c(c(1:22),"X","Y")
	gtf <- gtf[gtf$V1 %in% chr_index,]
	
	df_immu_region <- gtf[gtf$V2 %in% immu_id,]
	bed_immu_region <- data.frame(df_immu_region$V1,df_immu_region$V4, df_immu_region$V5)
	names(bed_immu_region) <- c("chr","start","end")
	bed_immu_region$chr <- paste("chr",as.character(bed_immu_region$chr),sep = "")
	write.table(bed_immu_region,"immu_region_from_gtf.bed", sep = "\t", row.names = F, col.names = F,quote = F)
	
	###将免疫高变区merge
	system("sort -k 1,1 -k 2,2n immu_region_from_gtf.bed > sort_chr_immu_region_from_gtf.bed")
	system("bedtools merge -i sort_chr_immu_region_from_gtf.bed > merged_sort_chr_immu_region_from_gtf.bed")

	
	###将免疫高变区前后扩增ext KB
	ext = 50000
	merged <- fread("merged_sort_chr_immu_region_from_gtf.bed")
	merged$V2 <- merged$V2 - ext
	merged$V3 <- merged$V3 + ext
	write.table(merged,"extended_merged_sort_chr_immu_region_from_gtf.bed", sep = "\t", row.names = F ,col.names = F, quote = F)

	###再次merge
	system("sort -k 1,1 -k 2,2n extended_merged_sort_chr_immu_region_from_gtf.bed > sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
	system("bedtools merge -i sorted_extended_merged_sort_chr_immu_region_from_gtf.bed > merged_sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
	
	file.remove("sort_chr_immu_region_from_gtf.bed",
				"merged_sort_chr_immu_region_from_gtf.bed",
				"immu_region_from_gtf.bed",
				"sorted_extended_merged_sort_chr_immu_region_from_gtf.bed",
				"extended_merged_sort_chr_immu_region_from_gtf.bed")

}				


###获取免疫和外显子两个区域的并集区间
{
	cds <- fread("cds_region.bed")
	immu <- fread("merged_sorted_extended_merged_sort_chr_immu_region_from_gtf.bed")
	result <- rbind(cds,immu)
	write.table(result, "rbind_three.bed", sep = "\t", row.names = F, col.names = F, quote = F)	
	
	system("sort -k 1,1 -k 2,2n rbind_three.bed > sorted_rbind_three.bed")
	system("bedtools merge -i sorted_rbind_three.bed > merged_sorted_rbind_three.bed")
	
	file.remove("rbind_three.bed","sorted_rbind_three.bed")
}


###获取基因组上除去免疫和外显子的剩下的非编码区间
{
	system("sort -k 1,1 -k 2,2n merged_sorted_rbind_three.bed > sorted_merged_sorted_rbind_three.bed")
	system("sort -k 1,1 -k 2,2n human.hg19.genome > sorted_human.hg19.genome")
	######获取全基因组mask区间的补集######
	system("bedtools complement -i sorted_merged_sorted_rbind_three.bed -g sorted_human.hg19.genome > hg19_nonmask_region.bed")
	
	file.remove("sorted_merged_sorted_rbind_three.bed","sorted_human.hg19.genome")
}



###获取非编码区的所有体细胞突变
{
options(scipen = 30)
df <- fread("rmSNP_single_mut.bed")
write.table(df,"rmSNP_single_mut.bed", sep = "\t", row.names = F, quote = F, col.names = F)

region_mut("rmSNP_single_mut.bed","hg19_nonmask_region.bed","noncoding_mut")
}









	





	
	
	







