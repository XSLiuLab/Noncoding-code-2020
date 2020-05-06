library(dplyr)

# loc = dplyr::tibble(
#   chrom = "chr1",
#   start = 43824626 - 3000 - 1,
#   end = 43824626 + 3000
# )

## Use mutaiton motif as center
loc = dplyr::tibble(
  chrom = "chr1",
  start = 43824527 - 500,
  end = 43824527 + 500
)

# Signal from ENCODE------------------------------------------------------------------

HEK293 = rtracklayer::import.bw("~/biodata/CDC20/ENCFF000WXM.bigWig")
HEK293 = data.table::as.data.table(HEK293)
HEK293_CDC20 = HEK293[start < loc$end & end > loc$start & seqnames == loc$chrom]
summary(HEK293_CDC20$score)

HELA = rtracklayer::import.bw("~/biodata/CDC20/ENCFF000XDZ.bigWig")
HELA = data.table::as.data.table(HELA)
HELA_CDC20 = HELA[start < loc$end & end > loc$start & seqnames == loc$chrom]
summary(HELA_CDC20$score)

save(HEK293_CDC20, HELA_CDC20, file = "signal_plotting.RData")

# Plotting ----------------------------------------------------------------
load(file = "signal_plotting.RData")
library(ggplot2)
library(ggforce)

start_point = 43824626
motif_start = 43824524
motif_end = 43824529

dat = rbind(
  HEK293_CDC20[, c("start", "score")],
  HEK293_CDC20[, c("end", "score")],
  use.names = FALSE
)
dat = unique(dat)

p1 <- ggplot(dat, aes(x=start, y=score)) + 
  geom_area() +
  geom_vline(xintercept = start_point, color = "green") + 
  annotate("rect", xmin = motif_start, xmax = motif_end, ymin = 0, ymax = 15, alpha = .3, fill = "red") + 
  cowplot::theme_cowplot() + 
  labs(x = NULL, y = "Signal in HEK293 cell line")

## Add zoom in
# ggplot(dat, aes(x=start, y=score)) + 
#   geom_area() +
#   # xlim(start_point - 100, start_point + 10) +
#   geom_vline(xintercept = start_point, color = "green") + 
#   facet_zoom(xlim = c(start_point - 150, start_point + 100)) +
#   cowplot::theme_cowplot() + 
#   labs(x = "Genome coordinates of chr1", y = "Chip-Seq signal score in HEK293 cell line")

# dat[order(start)][start > start_point + 120]

dat2 = rbind(
  HELA_CDC20[, c("start", "score")],
  HELA_CDC20[, c("end", "score")],
  use.names = FALSE
)
dat2 = unique(dat2)

p2 <- ggplot(dat2, aes(x=start, y=score)) +
  geom_area() +
  geom_vline(xintercept = start_point, color = "green") +
  annotate("rect", xmin = motif_start, xmax = motif_end, ymin = 0, ymax = 8.4, alpha = .3, fill = "red") + 
  cowplot::theme_cowplot() +
  labs(x = "Genome coordinates of chr1", y = "Signal score in HELA cell line")

library(patchwork)
p = p1 / p2

ggsave(filename = "ELK4_coverage.pdf", plot = p, width = 6, height = 4)

## Add zoom in
# ggplot(dat2, aes(x=start, y=score)) +
#   geom_area() +
#   # xlim(start_point - 100, start_point + 10) +
#   geom_vline(xintercept = start_point, color = "green") +
#   facet_zoom(xlim = c(start_point - 150, start_point + 100)) +
#   cowplot::theme_cowplot() +
#   labs(x = "Genome coordinates of chr1", y = "Chip-Seq signal score in HELA cell line")

# Read coverage -----------------------------------------------------------

# readr:::write_tsv(loc, path = "cdc20_selected_region.bed", col_names = FALSE)

# 
# getCoverage = function(bam_file, bed_file) {
#   stopifnot(length(bam_file) == 1)
#   out_file = basename(sub("bam", "coverage", bam_file))
#   # system(paste(
#   #   "/Users/wsx/anaconda3/envs/biosoft/bin/samtools depth",
#   #   "-b",
#   #   bed_file,
#   #   bam_file,
#   #   ">",
#   #   out_file
#   # ))
#   
#   # ref: https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
#   system(
#     paste(
#       "/Users/wsx/anaconda3/envs/biosoft/bin/bedtools coverage",
#       "-a", bed_file,
#       "-b", bam_file,
#       "-d",
#       ">", out_file
#     )
#   )
#   
# }
# 
# getCoverage(bam_file = "~/biodata/CDC20/ENCFF000WXK.bam", bed_file = "cdc20_selected_region.bed")
# getCoverage(bam_file = "~/biodata/CDC20/ENCFF000WXN.bam", bed_file = "cdc20_selected_region.bed")
# 
# dt <- readr::read_tsv("ENCFF000WXN.coverage", col_names = F)


