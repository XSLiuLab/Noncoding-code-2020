library(tidyverse)
library(data.table)
library(ggrepel)

load(file = "RegionMutationList.RData")

prob_region_final[, p_value := -log10(p_val)][, p_value := ifelse(p_value > 15, 15, p_value)]

pos <- position_jitter(width = 0.2, height = 0.03, seed = 2)
p1 = ggplot(prob_region_final[gene_name != "BCL2"], aes(x = count, y = p_value)) + 
  geom_jitter(color = ifelse(prob_region_final[gene_name != "BCL2"]$count > 20, "red", "black"),
              position = pos, size = 1, alpha = 0.7) +
  geom_text_repel(data = dplyr::filter(prob_region_final, count > 20 & gene_name != "BCL2"),
                   aes(label = gene_name),
                  hjust = 1,
                  position = pos,
                  segment.size = 0.4,min.segment.length = 0.1,
                  size = 2) +
  labs(x = "Number of samples with base mutated", y = "-log10 P value") +
  cowplot::theme_cowplot(rel_small = 1.5) 

p1

ggsave(filename = "pvalue_vs_count.pdf", plot = p1, width = 6, height = 5)


# p1 = ggplot(prob_region_final, aes(x = count, y = p_value)) + 
#   geom_jitter(color = ifelse((prob_region_final$count > 10 & prob_region_final$p_value > 5), "red", 
#                              ifelse(prob_region_final$count >= 30, "blue", "black"))) +
#   geom_label_repel(data = dplyr::filter(prob_region_final, (count > 10 & p_value > 5) | count >= 30 ),
#                    aes(label = gene_name),
#                    angle        = 0,
#                    #vjust        = 1,
#                    #hjust = 1,
#                    segment.size = 0.2) +
#   labs(x = "Number of samples with base mutated", y = "-log10 P value") +
#   cowplot::theme_cowplot(rel_small = 1.5) 