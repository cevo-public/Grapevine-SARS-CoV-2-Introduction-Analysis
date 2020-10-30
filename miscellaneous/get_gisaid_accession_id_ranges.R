# Get ranges of GISAID EPI ISLS from our seqs

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
MASTER_METADATA <- "/Volumes/nadeaus/2019-ncov-data/gisaid_data/2020-10-07_ch_cluster_analysis/nextmeta_with_unreleased.tsv"
OUTDIR <- paste(WORKDIR, "accession_ids", sep = "/")

system(command = paste("mkdir -p", OUTDIR))

master_metadata <- read.delim(file = MASTER_METADATA)

our_seq_metadata <- master_metadata %>% 
  filter(originating_lab == "Viollier AG")

accession_id_info <- data.frame(id = our_seq_metadata$gisaid_epi_isl) %>%
  tidyr::separate(col = id, into = c("p1", "p2", "number"), remove = F) %>%
  mutate(number = as.numeric(number)) %>%
  arrange(number) %>%
  filter(id != "TBD")

is_first <- T
starts <- c()
ends <- c()
prev_number <- accession_id_info[1, "number"]
for (i in 1:nrow(accession_id_info)) {
  current_number <- accession_id_info[i, "number"]
  if (is_first) {
    is_first <- F
    starts <- c(starts, paste("EPI_ISL_", prev_number, sep = ""))
  } else if (current_number != prev_number + 1) {
    ends <- c(ends, paste("EPI_ISL_", prev_number, sep = ""))
    is_first <- T
  }
  prev_number <- current_number
}  
ends <- c(ends, paste("EPI_ISL_", current_number, sep = ""))

accession_id_ranges <- data.frame(starts, ends)

write.table(
  file = paste(OUTDIR, "gisaid_accession_id_ranges.txt", sep = "/"),
  x = accession_id_ranges, row.names = F, col.names = T, quote = F, sep = "\t")
