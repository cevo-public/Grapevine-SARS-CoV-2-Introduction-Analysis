# This script identifies and formats sequences that have been QC'd but aren't 
# included in the nextfasta download from GISAID yet so they can be concatenated
# with the nextfasta data.

# Note that ETH ID 240000 never got released (I think because it wasn't in the metadata?) so it will be tacked on with any not-yet included seqs.

# OUR_SEQS <- "/Volumes/nadeaus/2019-ncov-data/gisaid_data/2020-10-07_ch_cluster_analysis/our_qcd_seqs.fasta"
# OUR_SEQ_METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/ViollierConfidential/Metadata/merged_metadata/merged_metadata_14.txt"
# NEXTMETA <- "/Volumes/nadeaus/2019-ncov-data/gisaid_data/2020-10-07_ch_cluster_analysis/nextmeta.tsv"
# OUTDIR <- "~/Downloads"
# UTILITY_FUNCTIONS <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis/scripts/utility_functions.R"

require(ape)
require(argparse)
require(dplyr)

parser <- argparse::ArgumentParser()
parser$add_argument("--ourseqsreleased", type="character")
parser$add_argument("--ourseqsmetadata", type="character")
parser$add_argument("--nextmeta", type="character")
parser$add_argument("--outdir", type="character")
parser$add_argument("--utilityfunctions", type="character")

args <- parser$parse_args()

OUR_SEQS <- args$ourseqsreleased
OUR_SEQ_METADATA <- args$ourseqsmetadata
NEXTMETA <- args$nextmeta
OUTDIR <- args$outdir
UTILITY_FUNCTIONS <- args$utilityfunctions

source(UTILITY_FUNCTIONS)
system(command = paste("mkdir -p", OUTDIR))

# Load data
nextmeta <- read.delim(file = NEXTMETA, sep = "\t", stringsAsFactors = F)
our_metadata <- read.delim(file = OUR_SEQ_METADATA, sep = ";", stringsAsFactors = F)
our_seqs <- ape::read.FASTA(file = OUR_SEQS)

# Find samples that haven't been released yet but pass QC ----------------------

our_seqs_nextmeta <- nextmeta %>%
  filter(authors == "Christian Beisel et al")
get_eth_id_from_strain_name <- function(strain_name) {
  removed_prefix <- strsplit(strain_name, split = "-")[[1]][3]
  removed_suffix <- strsplit(removed_prefix, split = "/")[[1]][1]
  removed_tail <- strsplit(removed_suffix, split = "_")[[1]][1]
  removed_plus <- strsplit(removed_tail, split = "plus")[[1]][1]
  return(removed_plus)
}

nextmeta_eth_ids <- unlist(lapply(X = our_seqs_nextmeta$strain, FUN = get_eth_id_from_strain_name))

our_seq_names <- names(our_seqs)
get_eth_id_from_sample_name <- function(strain_name) {
  if(grepl(pattern = "^[[:alpha:]]", x = strain_name)) {  # ZRH sample name format
    return(strsplit(x = strain_name, split = "_")[[1]][2])
  }
  removed_tail <- strsplit(strain_name, split = "_")[[1]][1]
  return(removed_tail)
}
our_seqs_eth_ids <- unlist(lapply(X = our_seq_names, FUN = get_eth_id_from_sample_name))

eth_ids_not_in_nextmeta <- setdiff(our_seqs_eth_ids, nextmeta_eth_ids)
eth_ids_not_in_qc_dir <- setdiff(nextmeta_eth_ids, our_seqs_eth_ids)

if (length(eth_ids_not_in_qc_dir) > 0) {
  print(eth_ids_not_in_qc_dir)
  warning(paste(length(eth_ids_not_in_qc_dir), "ETH IDs identified from Nextdata download not in ETH IDs identified from our QC'd sequences."))
}

# Write out fasta file of our sequences to include -----------------------------
names(our_seqs) <- our_seqs_eth_ids
our_seqs_to_add <- subset(
  x = our_seqs,
  subset = names(our_seqs) %in% eth_ids_not_in_nextmeta)

# Need to add Switzerland for grepping of focal Swiss seqs later
eth_ids_added <- names(our_seqs_to_add)
names(our_seqs_to_add) <- paste("Switzerland_XX-ETHZ-", names(our_seqs_to_add), "_2020", sep = "")

print(paste("Writing out", length(our_seqs_to_add), "sequences to add."))

ape::write.FASTA(
  x = our_seqs_to_add,
  file = paste(OUTDIR, "our_seqs_to_add_to_nextdata.fasta", sep = "/"))

# Write out metadata file of our sequences to include --------------------------

our_seqs_to_add_metadata <- our_metadata %>%
  filter(ETH.ID %in% eth_ids_added) %>%
  mutate(
    division = unlist(lapply(X = Kanton, FUN = get_ns_division_from_canton_code)))
  
seq_lengths <- unlist(lapply(X = our_seqs_to_add, FUN = length))
if (!all(seq_lengths == 29903)) {
  stop("Not all sequences have expected length 29903.")
} 

our_seqs_to_add_metadata_ns_format <- data.frame(
  strain = paste("Switzerland_XX-ETHZ-", our_seqs_to_add_metadata$ETH.ID, "_2020", sep = ""),
  virus = "ncov",
  gisaid_epi_isl = "TBD",
  genbank_accession = "?",
  date = as.Date(our_seqs_to_add_metadata$Order.date, format = "%Y-%m-%d"),
  region = "Europe",
  country = "Switzerland",
  division = our_seqs_to_add_metadata$division,
  location = "",
  region_exposure = "Europe",
  country_exposure = "Switzerland",
  division_exposure = our_seqs_to_add_metadata$division,
  segment = "genome",
  length = 29903,
  host = "Human",
  age = "?",
  sex = "?",
  pangolin_lineage = "TBD",
  GISAID_clade = "TBD",
  originating_lab = "Viollier AG",
  submitting_lab = "Department of Biosystems Science and Engineering, ETH Zürich",
  authors = "Christian Beisel et al",
  url = "?",
  title = "?",
  paper_url = "?",
  date_submitted = "TBD")

write.table(
  x = our_seqs_to_add_metadata_ns_format, 
  file = paste(OUTDIR, "our_seqs_to_add_to_nextdata.tsv", sep = "/"),
  sep = "\t", col.names = F, row.names = F, quote = F)







