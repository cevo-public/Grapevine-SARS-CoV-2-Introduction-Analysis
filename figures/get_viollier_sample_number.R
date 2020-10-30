OUTDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/figures"

master_metadata <- read.delim(
  file = "/Volumes/nadeaus/2019-ncov-data/gisaid_data/2020-10-07_ch_cluster_analysis/nextmeta_with_unreleased.tsv")
tree_metadata_1 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_1_n_sim_1000_n_imports_padded_0/rep_1_n_sim_1000_n_imports_padded_0_tree_metadata.txt",
  quote = "")
tree_metadata_2 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_2_n_sim_1000_n_imports_padded_0/rep_2_n_sim_1000_n_imports_padded_0_tree_metadata.txt",
  quote = "")
tree_metadata_3 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_3_n_sim_1000_n_imports_padded_0/rep_3_n_sim_1000_n_imports_padded_0_tree_metadata.txt",
  quote = "")
context_metadata_1 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_1_n_sim_1000_n_imports_padded_0/rep_1_n_sim_1000_n_imports_padded_0_context_metadata.txt",
  quote = "")
context_metadata_2 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_2_n_sim_1000_n_imports_padded_0/rep_2_n_sim_1000_n_imports_padded_0_context_metadata.txt",
  quote = "")
context_metadata_3 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_3_n_sim_1000_n_imports_padded_0/rep_3_n_sim_1000_n_imports_padded_0_context_metadata.txt",
  quote = "")
similarity_metadata_1 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_1_n_sim_1000_n_imports_padded_0/rep_1_n_sim_1000_n_imports_padded_0_priority_metadata.txt",
  quote = "")
similarity_metadata_2 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_2_n_sim_1000_n_imports_padded_0/rep_2_n_sim_1000_n_imports_padded_0_priority_metadata.txt",
  quote = "")
similarity_metadata_3 <- read.delim(
  file = "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/data/alignments/rep_3_n_sim_1000_n_imports_padded_0/rep_3_n_sim_1000_n_imports_padded_0_priority_metadata.txt",
  quote = "")

mm_v <- master_metadata %>% filter(originating_lab == "Viollier AG", as.Date(date) < as.Date("2020-08-31"))  # Augur filter is end exclusive, so seqs ON 31. Aug not included
nrow(tree_metadata_1 %>% filter(originating_lab == "Viollier AG"))
nrow(tree_metadata_2 %>% filter(originating_lab == "Viollier AG"))
nrow(tree_metadata_2 %>% filter(originating_lab == "Viollier AG"))

nrow(tree_metadata_1 %>% filter(country_recoded == "Switzerland", originating_lab != "Viollier AG"))
nrow(tree_metadata_2 %>% filter(country_recoded == "Switzerland", originating_lab != "Viollier AG"))
nrow( tree_metadata_2 %>% filter(country_recoded == "Switzerland", originating_lab != "Viollier AG"))

setdiff(mm_v$strain, tm1_v$old_strain)
setdiff(mm_v$strain, tm2_v$old_strain)
setdiff(mm_v$strain, tm3_v$old_strain)  # these should all be the same

# Final dataset numbers
nrow(context_metadata_1)
nrow(context_metadata_2)
nrow(context_metadata_3)

setdiff(context_metadata_1$strain, context_metadata_2$strain)
setdiff(context_metadata_1$strain, context_metadata_3$strain)
setdiff(context_metadata_2$strain, context_metadata_3$strain)  # these should all be non-empty

nrow(similarity_metadata_1)
nrow(similarity_metadata_2)
nrow(similarity_metadata_3)

setdiff(similarity_metadata_1$strain, similarity_metadata_2$strain)
setdiff(similarity_metadata_1$strain, similarity_metadata_3$strain)
setdiff(similarity_metadata_2$strain, similarity_metadata_3$strain)  # these should all be empty

nrow(tree_metadata_1)
nrow(tree_metadata_2)
nrow(tree_metadata_3)  # slightly diff b/c diff overlap b/w similarity and context sequences

# Make GISAID acknowledgments table
seqs_analyzed <- master_metadata %>% filter(
  gisaid_epi_isl %in% tree_metadata_1$gisaid_epi_isl | 
    gisaid_epi_isl %in% tree_metadata_2$gisaid_epi_isl | 
    gisaid_epi_isl %in% tree_metadata_3$gisaid_epi_isl)

gisaid_table <- seqs_analyzed %>% select(
  gisaid_epi_isl, originating_lab, submitting_lab, authors) %>%
  group_by(originating_lab, submitting_lab, authors) %>%
  summarize("Accession IDs" = paste0(gisaid_epi_isl, collapse = ", ")) %>%
  tidyr::unite("details", originating_lab, submitting_lab, authors, sep = "\t") %>%
  select(details, `Accession IDs`)

write.table(
  x = as.vector(t(as.matrix(gisaid_table))), 
  file = paste(OUTDIR, "gisaid_table.txt", sep = "/"),
  row.names = F, quote = F, sep = "\t")

