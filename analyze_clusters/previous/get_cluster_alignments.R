# Output fasta file with unique cluster identifier in header.

require(ape)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_3_n_sim_1000_n_imports_padded_0"
PREFIX <- "rep_3_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F"

CLUSTERS <- paste(WORKDIR, "/clusters/", PREFIX, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
OUTDIR <- paste(WORKDIR, "/data/cluster_alignments/", sep = "")
SWISS_ALIGNMENT <- paste(WORKDIR, "/data/qc_master_alignment/swiss_alignment_filtered2_masked_oneline.fasta", sep = "")

source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

# Load data
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F) %>% unique.data.frame()  # remove duplicates (to account for a concatenate mistake in metadata)
cluster_data <- read.delim(file = CLUSTERS)
swiss_alignment <- ape::read.FASTA(file = SWISS_ALIGNMENT)

cluster_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = cluster_data, metadata = metadata) %>%
  filter(originating_lab == "Viollier AG")

# Generate new sequence headers of form:
# <strain with "\" & " " --> "_">|<gisaid_epi_isl>|<yyyy-mm-dd date>|<cluster_idx>
# View(names(swiss_alignment))

ch_seq_metadata <- data.frame(old_strain = names(swiss_alignment))
ch_seq_metadata_2 <- merge(
  x = ch_seq_metadata, y = metadata[c("old_strain", "strain")],
  all.x = T)

ch_seq_metadata_final <- merge(
  x = ch_seq_metadata_2, y = cluster_data_by_tip[c("tip", "cluster_idx")],
  by.x = "strain", by.y = "tip",
  all.y = T)  # some seqs thrown out be least-squares-dating for violating strict clock model

# ggplot(
#   data = ch_seq_metadata_final,
#   aes(x = cluster_idx)) +
#   geom_bar()

ch_seq_metadata_final$new_header <- unlist(mapply(
  FUN = paste,
  ch_seq_metadata_final$strain,
  ch_seq_metadata_final$cluster_idx, 
  sep = "|"))

# Apply new sequence headers
swiss_alignment_cluster_seqs <- subset(
  x = swiss_alignment, 
  subset = names(swiss_alignment) %in% ch_seq_metadata_final$old_strain)

table_idxs <- match(names(swiss_alignment_cluster_seqs), ch_seq_metadata_final$old_strain)
new_headers <- ch_seq_metadata_final[table_idxs, "new_header"]
cat(paste0(head(new_headers), collapse = "\n"))
cat(paste(head(names(swiss_alignment_cluster_seqs))), sep = "\n")
cat(paste0(tail(new_headers), collapse = "\n"))
cat(paste(tail(names(swiss_alignment_cluster_seqs))), sep = "\n")

names(swiss_alignment_cluster_seqs) <- new_headers
ape::write.FASTA(
  x = swiss_alignment_cluster_seqs,
  file = paste(OUTDIR, "/", PREFIX, "_viollier_only_swiss_clusters.fasta", sep = ""))
