# Generate table for text with cluster stats from different analyses

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
DATADIR <- paste(WORKDIR, "cluster_stats", sep = "/")
OUTDIR <- paste(WORKDIR, "/figures/cluster_stats", sep = "")

system(command = paste("mkdir -p", OUTDIR))

stat_files <- list.files(
  path = DATADIR, full.names = F, pattern = "_cluster_stats.txt")

is_first <- T
for (file in stat_files) {
  stats <- read.delim(file = paste(DATADIR, file, sep = "/"))
  prefix <- gsub(x = file, pattern = "_cluster_stats.txt", replacement = "")
  M <- strsplit(x = prefix, split = "_")[[1]][11]
  P <- strsplit(x = prefix, split = "_")[[1]][13]
  polyswiss <- strsplit(x = prefix, split = "_")[[1]][16]
  stats$analysis <- paste0(x = strsplit(prefix, split = "_")[[1]][1:2], collapse = "_")
  stats$M <- M
  stats$P <- P
  stats$polyswiss <- polyswiss
  if (is_first) {
    all_stats <- stats
    is_first <- F
  } else {
    all_stats <- rbind(all_stats, stats)
  }
}

all_stats <- all_stats %>% arrange(P)

pos1 <- position_jitter(width = 0.15, height = 0, seed = 1)
p1 <- ggplot(
  data = all_stats,
  aes(x = P, y = n_transmission_chain, label = M)) +
  geom_point(position = pos1) + 
  geom_text_repel(position = pos1) +
  # geom_text(position = position_jitter(width = 0.15)) + 
  labs(x = element_blank(), y = "No. transmission\nchains") + 
  theme_bw() +
  theme(legend.position = "none")
show(p1)

pos2 <- position_jitter(width = 0.15, height = 0, seed = 1)
p2 <- ggplot(
  data = all_stats,
  aes(x = P, y = n_singletons, label = M)) +
  geom_point(position = pos2) + 
  geom_text_repel(position = pos2) +
  # geom_text(position = position_jitter(width = 0.15)) + 
  labs(x = element_blank(), y = "No. singletons") + 
  theme_bw() +
  theme(legend.position = "none")

pos3 <- position_jitter(width = 0.15, height = 0, seed = 1)
p3 <- ggplot(
  data = all_stats,
  aes(x = P, y = biggest_cluster, label = M)) +
  geom_point(position = pos3) + 
  geom_text_repel(position = pos3) +
  # geom_text(position = pos, aes(vjust = 0.5)) +
  labs(x = element_blank(), y = "Largest cluster") + 
  theme_bw() + 
  theme(legend.position = "bottom")
show(p3)

library(gtable)
library(grid)
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)

png(
  file = paste(OUTDIR, "cluster_stats.png", sep = "/"),
  width = 6.5, height = 4.25, units = "in", res = 300)
grid.newpage()
grid.draw(g)
dev.off()

write.table(
  x = all_stats,
  file = paste(OUTDIR, "cluster_stats.txt", sep = "/"), 
  quote = F, row.names = F, col.names = T, sep = "\t")
