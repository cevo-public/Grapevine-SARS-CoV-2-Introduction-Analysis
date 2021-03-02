# Plot # seqs in the tree from each country vs. source location contribution to 
# see if any countries stand out in particular.

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-29_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
PREFIX <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F"

METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
CONTEXT_METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_context_metadata.txt", sep = "")
TREE_DATA_WITH_ASR <- paste(WORKDIR, "/asr/", PREFIX, "_tree_data_with_asr.txt", sep = "")
CLUSTER_DATA <- CLUSTER_DATA <- paste(WORKDIR, "/clusters/", PREFIX, "_clusters.txt", sep = "")
OUTDIR <- paste(WORKDIR, "figures/asr_contributions", sep = "/")

source(paste(WORKDIR, "scripts/utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

# Load data
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
context_metadata <- read.table(file = CONTEXT_METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
tree_data_with_asr <- read.delim(file = TREE_DATA_WITH_ASR, stringsAsFactors = F)
cluster_data <- read.delim(file = CLUSTER_DATA, stringsAsFactors = F)

# Get No. samples per country per month
context_samples_per_month <- context_metadata %>%
  mutate(month = lubridate::month(date)) %>%
  group_by(country_recoded, month) %>%
  summarize(n_seqs = n()) %>%
  rename(country = country_recoded)

# Get ASR contribution per country per month
cluster_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = cluster_data, metadata = metadata)

mrca_asr_data_all <- cluster_data_by_tip %>%
  group_by(foreign_mrca, foreign_tmrca) %>%
  summarise(n_swiss_descendents = sum(size))

mrca_asr_data_all <- merge(
  x = mrca_asr_data_all, y = tree_data_with_asr,
  by.x = "foreign_mrca", by.y = "node", all.x = T)

loc_colnames <- colnames(mrca_asr_data_all)[grep(
  x = colnames(mrca_asr_data_all),
  pattern = "_loc_weight")]
loc_colnames <- unlist(lapply(X = loc_colnames, FUN = recode_colnames))

# Recode ASR locations (format "First.Second_loc_weight") to match country_recoded names
colnames(mrca_asr_data_all) <- unlist(lapply(
  X =  colnames(mrca_asr_data_all), 
  FUN = recode_colnames))

# Make long data frame for summary
mrca_asr_data_all_long <- mrca_asr_data_all %>% tidyr::pivot_longer(
  cols = all_of(loc_colnames),
  names_to = "country",
  values_to = "asr_contribution") %>% 
  mutate(
    asr_contribution = as.numeric(as.character(asr_contribution)),
    month = lubridate::month(as.Date(foreign_tmrca)))

asr_contribution_by_month <- mrca_asr_data_all_long %>%
  group_by(country, month) %>%
  summarize(asr_contribution = sum(asr_contribution, na.rm = T))

asr_contribution_vs_seqs_per_month <- merge(
  x = context_samples_per_month, y = asr_contribution_by_month,
  all = T)
asr_contribution_vs_seqs_per_month[is.na(asr_contribution_vs_seqs_per_month$asr_contribution), "asr_contribution"] <- 0
asr_contribution_vs_seqs_per_month[is.na(asr_contribution_vs_seqs_per_month$n_seqs), "n_seqs"] <- 0

asr_contribution_vs_seqs_per_month <- asr_contribution_vs_seqs_per_month %>%
  filter(!(country %in% c("dummy_loc", "Switzerland")))

asr_contribution_vs_seqs_per_month <- asr_contribution_vs_seqs_per_month %>% mutate(
  label = case_when(
    asr_contribution > 1 | n_seqs > 2 ~ paste(country, month, sep = ": "),
    T ~ ""))

asr_contribution_vs_seqs_per_month_plot <- ggplot(
  data = asr_contribution_vs_seqs_per_month,
  aes(x = n_seqs, y = asr_contribution)) + 
  geom_point() + 
  geom_text_repel(aes(label = label), size = 1.5) +
  scale_x_continuous(
    trans=scales::pseudo_log_trans(base = 10)) + 
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10)) + 
  theme_bw() + 
  labs(x = "No. context sequences", y = "Source location contribution")

png(
  file = paste(OUTDIR, paste(PREFIX, "asr_contribution_vs_seqs_per_month.png", sep = "_"), sep = "/"), 
  width = 4, height = 3, res = 300, units = "in")
print(asr_contribution_vs_seqs_per_month_plot)
dev.off()

asr_contribution_vs_total_seqs <- asr_contribution_vs_seqs_per_month %>%
  group_by(country) %>%
  summarize(asr_contribution_total = sum(asr_contribution), n_seqs_total = sum(n_seqs))

ggplot(
  data = asr_contribution_vs_total_seqs,
  aes(x = n_seqs_total, y = asr_contribution_total)) + 
  geom_point() + 
  geom_text_repel(aes(label = country), size = 2.5) +
  scale_x_continuous(
    trans=scales::pseudo_log_trans(base = 10)) + 
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10)) + 
  theme_bw() + 
  labs(x = "No. context sequences", y = "Source location contribution")

asr_contribution_vs_total_seqs_merged <- merge(
  x = asr_contribution_vs_total_seqs, y = asr_contribution_vs_seqs_per_month)

ggplot(
  data = asr_contribution_vs_total_seqs_merged %>% filter(month == 9),
  aes(x = n_seqs_total, y = asr_contribution)) + 
  geom_point() + 
  geom_text_repel(aes(label = country), size = 2.5) +
  scale_color_distiller(palette = "Spectral") +  
  scale_x_continuous(
    trans=scales::pseudo_log_trans(base = 10)) + 
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10)) + 
  theme_bw() + 
  labs(x = "Total context sequences in tree", y = "Source location contribution in month")

png(
  file = paste(OUTDIR, paste(PREFIX, "asr_contribution_vs_seqs_per_month.png", sep = "_"), sep = "/"), 
  width = 4, height = 3, res = 300, units = "in")
print(asr_contribution_vs_seqs_per_month_plot)
dev.off()

