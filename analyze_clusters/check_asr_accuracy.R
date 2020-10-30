# Check whether exposure information supports our cluster ASRs

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
PREFIX <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F"

TREE_DATA_WITH_ASR <- paste(WORKDIR, "/asr/", PREFIX, "_tree_data_with_asr.txt", sep = "")
CLUSTERS <- paste(WORKDIR, "/clusters/", PREFIX, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
OUTDIR <- paste(WORKDIR, "/figures/check_asr/", sep = "")

source(paste(WORKDIR, "scripts/utility_functions.R", sep = "/"))
source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

# Load data
metadata <- read.table(
  file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
cluster_data <- read.delim(file = CLUSTERS)
tree_data_with_asr <- read.delim(file = TREE_DATA_WITH_ASR, stringsAsFactors = F)

# Get cluster information ------------------------------------------------------
cluster_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = cluster_data, metadata = metadata)

cluster_data_by_tip_with_asr <- merge(
  x = cluster_data_by_tip, y = tree_data_with_asr,
  by.x = "foreign_mrca", by.y = "node")

# summarize_exposure_locs

cluster_data_by_tips_w_exposure <- cluster_data_by_tip_with_asr %>%
  group_by(foreign_mrca, cluster_idx) %>%
  filter(any(country != country_exposure)) 

cluster_data_by_tips_w_exposure_2 <- cluster_data_by_tips_w_exposure %>% pivot_longer(
  cols = all_of(asr_loc_colnames),
  names_to = "asr_loc", values_to = "asr_contribution")

summarize_exposure_locs <- function(country_exposure) {
  tab <- table(country_exposure)
  country_names <- names(tab)
  country_values <- unlist(tab)
  is_first <- T
  for (country in country_names) {
    if (is_first) {
      str <- paste(country, country_values[[country]], sep = ": ")
      is_first <- F
    } else {
      str <- paste(
        str,
        paste(country, country_values[[country]], sep = ": "),
        sep = ", ")
    }
  }
  return(str)
}

exposures_per_cluster <- cluster_data_by_tips_w_exposure %>% 
  group_by(cluster_idx) %>%
  summarize(
    n_unique_exposure_locs = length(unique(country_exposure[country_exposure != "Switzerland"])),
    n_exposed_samples = sum(country_exposure != country),
    exposure_info = summarize_exposure_locs(country_exposure))

write.table(
  x = exposures_per_cluster, 
  file = paste(OUTDIR, "/", PREFIX, "_exposures_per_cluster.txt", sep = ""),
  row.names = F, col.names = T, quote = F, sep = "\t")

asr_loc_colnames <- grep(x = colnames(tree_data_with_asr), value = T, pattern = "loc_weight")
summarize_asr_info <- function(node, tree_data_with_asr) {
  loc_info <- tree_data_with_asr[tree_data_with_asr$node == node, asr_loc_colnames]
  loc_info <- sort(loc_info, decreasing = T)
  nonzero_loc_values <- loc_info[loc_info > 0 & !is.na(loc_info)]
  nonzero_loc_names <- names(loc_info)[loc_info > 0 & !is.na(loc_info)]
  is_first <- T
  for (i in 1:length(nonzero_loc_names)) {
    if (is_first) {
      str <- paste(nonzero_loc_names[i], nonzero_loc_values[i], sep = ": ")
      is_first <- F
    } else {
      str <- paste(
        str,
        paste(nonzero_loc_names[i], nonzero_loc_values[i], sep = ": "),
        sep = ", ")
    }
  }
  return(str)
}

asr_inference_accuracy <- cluster_data_by_tips_w_exposure %>%
  group_by(foreign_mrca) %>%
  summarize(
    exposure_info = summarize_exposure_locs(country_exposure))

asr_inference_accuracy$asr_info <- unlist(lapply(
  X = asr_inference_accuracy$foreign_mrca, 
  FUN = summarize_asr_info,
  tree_data_with_asr = tree_data_with_asr))

write.table(
  x = asr_inference_accuracy, 
  file = paste(OUTDIR, "/", PREFIX, "_asr_accuracy.txt", sep = ""),
  row.names = F, col.names = T, quote = F, sep = "\t")