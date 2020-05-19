## Calculate p values for mutations and regions (freq>5)

# # Get the mutation frequency
# system(
#   "cat pred*.tsv.gz | grep -v donor | grep -v mock | cut -f2,3 | sort | uniq -c | sort -n | gzip > final_mutation_freq.tsv.gz"
# )
# 
# system(
#   "cat pred*.tsv.gz | grep -v donor | grep -v mock | gzip > final_mutation.tsv.gz"
# )

library(data.table)

freq_dt = fread("final_mutation_freq.tsv.gz", header = FALSE)
head(freq_dt)
freq_dt$chr = sub("[0-9]+ (chr[0-9XYxy]{1,2})", "\\1", freq_dt$V1)
table(freq_dt$chr)

freq_dt$freq = sub("([0-9]+) chr.*", "\\1", freq_dt$V1)
freq_dt$freq = as.integer(freq_dt$freq)

## Set the frequency threshold
freq_dt = freq_dt[freq > 3]
freq_dt$V1 = NULL
freq_dt$start = freq_dt$V2 - 5
freq_dt$end = freq_dt$V2 + 5
freq_dt$V2 = NULL
freq_dt = freq_dt[, c("chr", "start", "end", "freq"), with = F]
# must be unique
freq_dt$freq = NULL
freq_dt = unique(freq_dt)
saveRDS(freq_dt, file = "recurrent_region.rds")

## Load all mutations
mut_dt = fread("final_mutation.tsv.gz", header = FALSE)
colnames(mut_dt) = c("donor", "chr", "end", "prob")
# must be unique
mut_dt = unique(mut_dt)
mut_dt$start = mut_dt$end
mut_dt = mut_dt[, c("chr", "start", "end", "prob", "donor"), with=F]

str(head(mut_dt))

setkey(freq_dt, chr, start, end)
final_dt = foverlaps(mut_dt, freq_dt, type = "within")
final_dt = final_dt[!is.na(start)]

head(final_dt)
final_dt[, region_midpoint := paste(chr, as.integer(start + 5), sep = ":")]
final_dt$i.start <- NULL
final_dt[, mut_index := paste(chr, as.integer(i.end), sep = ":")]
final_dt$i.end = NULL

head(final_dt)

saveRDS(final_dt, file = "final_mutation_regions.rds")


library(data.table)
final_dt <- readRDS("final_mutation_regions.rds")

## Merge adjacent regions
library(IRanges)
region_dt = unique(final_dt[, .(chr, start, end, region_midpoint)])
merged_dt = region_dt[, data.table::as.data.table(reduce(IRanges(start, end))), by = .(chr)]

setkey(merged_dt, chr, start, end)
final_dt = foverlaps(final_dt, merged_dt, type = "within")
final_dt$i.start = NULL
final_dt$i.end = NULL
final_dt$region_midpoint = NULL
final_dt[, region_range := paste0(chr, ":", start, "-", end)]

library(magrittr)
table(final_dt$mut_index) %>% sort(decreasing = TRUE) %>% head()

## Add a weight to prob for each donor to get sample-specific prob as Prof.Liu and JingZhang devised.
#system("zcat final_mutation.tsv.gz | cut -f 1 | sort | uniq -c | sed 's/^[ \t]*//g' | sed 's/ /\t/g' > donor_noncoding_mut_freq.tsv")
donor_freq <- fread("donor_noncoding_mut_freq.tsv", header = FALSE)
colnames(donor_freq) <- c("freq", "donor")
donor_freq[, weight := freq / sum(freq)]

final_dt <- merge(final_dt, donor_freq, by = "donor", all.x = TRUE)
final_dt[, prob := prob * weight]
any(is.na(final_dt$prob))
final_dt[, c("freq", "weight") := NULL]

nrow(final_dt)
nrow(unique(final_dt))

final_dt = unique(final_dt)

## Get region prob
## CumSum(Pi) ~ 1-CumProd(1-Pi)

## Todo: prob should be probs in all samples (including non-mutated samples)
prob_region = final_dt[, .(prob = mean(prob) * width), 
         by = .(donor, region_range)][,
           .(p_val = 1 - poibin::ppoibin(.N - 1, prob),
             donor_list = paste(unique(donor), collapse = ","),
             count = .N),
           by = region_range
         ]

# barplot(dbinom(0:1000, 1000, 0.001), xlim = c(0, 10))

### Method based on point prob 
# region_df <- merge(unique(final_dt[, .(region_range, mut_index, donor)]),
#                    prob_point[, .(mut_index, p_val)], 
#                    by = "mut_index", all.x = TRUE)
# region_df
# 
# cal_region_p = function(p) {
#   1 - cumprod(1 - p)[length(p)]
# }
# 
# prob_region <- region_df[, .(p_val = cal_region_p(p_val),
#                              donor_list = paste(unique(donor), collapse = ","),
#                              mutation_list = paste(unique(mut_index), collapse = ",")),
#                          by = .(region_range)]
### Method based on point prob 

## Set 0 to minimal p value
prob_region[, p_val := ifelse(p_val < .Machine$double.xmin, 
                             .Machine$double.xmin,
                             p_val)]
prob_region = prob_region[order(p_val)]
prob_region


## Annotate records
# system(" cat ~/predict_prob/Homo_sapiens.GRCh37.75.gtf | grep gene | grep protein_coding > protein_coding_genes.tsv")
# gene_dt <- fread("/Volumes/pan/zhangjing/protein_coding_genes.tsv", header = F)
# gene_dt <- gene_dt[V3 == "gene"]
# 
# extract_col <- function(x, name) {
#   library(magrittr)
#   stringr::str_extract(x, paste0(name, " ([^;]+);")) %>%
#     stringr::str_remove(paste0(name, " ")) %>%
#     stringr::str_remove_all("\"") %>%
#     stringr::str_remove(";")
# }
# 
# gene_dt[, gene_name := extract_col(V9, "gene_name")]
# gene_df = gene_dt[, .(V1, V4, V5, V7, gene_name)]
# colnames(gene_df) = c("chr", "start", "end", "strand", "gene_name")
# save(gene_df, file = "gene_df.RData")

load(file = "gene_df.RData")

gene_df[, chr := paste0("chr", chr)]
gene_df[, gene_start := start]
gene_df[, gene_end := end]
gene_df[, `:=`(
  start = ifelse(strand == "+", gene_start - 5000, gene_end + 1), 
  end   = ifelse(strand == "+", gene_start - 1, gene_end + 5000)
  )]

prob_region = tidyr::separate(prob_region, col = "region_range", into = c("chr", "start", "end"))
prob_region = as.data.table(prob_region)
prob_region$start = as.integer(prob_region$start)
prob_region$end = as.integer(prob_region$end)

prob_region_final <- foverlaps(
  prob_region,
  gene_df,
  type = "any"
)

prob_region_final = prob_region_final[!is.na(gene_name)][
  , .(gene_name, chr, i.start, i.end, p_val, donor_list, count)][
    order(p_val)]
prob_region_final = prob_region_final[order(count, decreasing = TRUE)]
colnames(prob_region_final)[3:4] = c("start", "end")

save(prob_region_final, file = "RegionMutationList.RData")
load(file = "RegionMutationList.RData")

prob_region_final = prob_region_final[gene_name != "BCL2"]
prob_region_final[, p_value := -log10(p_val)]#[, p_value := ifelse(p_value > 15, 15, p_value)]
prob_region_final[, pos := paste0(chr, ":", start, "-", end)]

openxlsx::write.xlsx(prob_region_final[, list(pos, gene_name, count, p_value)], file = "RegionMutationList.xlsx")

