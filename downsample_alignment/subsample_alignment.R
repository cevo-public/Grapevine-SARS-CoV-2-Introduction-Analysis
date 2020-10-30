# Downsampling protocols:
# a) get list of most similar sequences
# b) get list of context sequences based on estimated # imports
# c) get 2 sequences that define the outgroup
# d) get all Swiss sequences
# c) down-sample alignment to the union of these sequence sets

require(tidyr)
require(dplyr)
require(ape)
require(ggplot2)
require(argparse)

# WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-15_ch_cluster_analysis"
# METADATA <- paste(WORKDIR, "data/metadata_all.txt", sep = "/")
# ALIGNMENT <- paste(WORKDIR, "data/alignment_filtered2_masked_oneline.fasta", sep = "/")
# N_MOST_SIMILAR_SEQS <- 1000
# N_CONTEXT_SEQ_DATA <- paste(WORKDIR, "data/samples_per_country_month.txt", sep = "/")
# PREFIX <- "test"
# OUTDIR <- "~/Downloads"

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character")
parser$add_argument("--alignment", type="character")
parser$add_argument("--nsimseqs", type="integer")
parser$add_argument("--ncontextsamples", type="character")
parser$add_argument("--prefix", type="character")
parser$add_argument("--outdir", type = "character")
parser$add_argument("--maxyearmonthdec", type="double")

args <- parser$parse_args()

METADATA <- args$metadata
ALIGNMENT <- args$alignment
N_MOST_SIMILAR_SEQS <- args$nsimseqs
N_CONTEXT_SEQ_DATA <- args$ncontextsamples
PREFIX <- args$prefix
OUTDIR <- args$outdir
MAX_YEAR_MONTH_DEC <- args$maxyearmonthdec

# Hardcoded things
OUTGROUP_IDS <- c("Wuhan/Hu-1/2019",  "Wuhan/WH01/2019")
seq_fn <- paste(PREFIX, "alignment.fasta", sep = "_")
seq_fp <- paste(OUTDIR, seq_fn, sep = "/")

# Functions
get_id_from_sample_name <- function(sample_name) {
  return(unlist(strsplit(sample_name, split = "\\|"))[2])
}

fetch_seqs_from_gisaid_alignment = function(gisaid_fp, seq_headers, new_headers, seq_fp) {
  con <- file(gisaid_fp, "r")
  con2 <- file(seq_fp, "w")
  n <- 0
  headers_found <- c()
  while (T) {
    lines = readLines(con, n = 2)
    n <- n + 1
    if (n %% 10000 == 0) {
      print(paste("Scanning fasta line:", n))
    }
    if (is.na(lines[1])) {
      break
    }
    if (lines[1] %in% seq_headers) {
      headers_found <- c(headers_found, lines[1])
      # Change sequence headers to be:
      # <strain with "\" & " " --> "_">|<gisaid_epi_isl>|<yyyy-mm-dd date>
      lines[1] <- new_headers[match(lines[1], seq_headers)]
      writeLines(text = lines, con = con2)
    }
  }
  close(con)
  close(con2)
  return(headers_found)
}

# Setup
if (dir.exists(OUTDIR)) {
  print(OUTDIR)
  stop("OUTDIR file with this name already exists.")
} else {
  system(paste("mkdir -p", OUTDIR))
  
  # Load data
  metadata <- read.delim(file = METADATA, stringsAsFactors = F, sep = "\t", quote = "")
  samples_per_country_month <- read.delim(file = N_CONTEXT_SEQ_DATA, stringsAsFactors = F, sep = "\t")
  
  metadata$year_month_dec <- as.numeric(gsub(x = metadata$year_month, pattern = "-", replacement = "."))
  metadata <- metadata %>% filter(year_month_dec <= MAX_YEAR_MONTH_DEC)
  
  # Select top priority sequences
  priority_seq_data <- metadata %>%
    arrange(desc(priority)) %>%
    mutate(ordering = 1:n()) %>%
    top_n(wt = ordering, n = -N_MOST_SIMILAR_SEQS)
  
  priority_seq_summary <- priority_seq_data %>%
    group_by(country_recoded, year_month) %>% 
    summarize(n_samples = n())
  
  write.table(
    x = priority_seq_data,
    file = paste(OUTDIR, paste(PREFIX, "priority_metadata.txt", sep = "_"), sep = "/"),
    row.names = F, col.names = T, quote = F, sep = "\t")
  
  # Select context sequences randomly based on specified # per country per month
  is_first <- T
  for (country_c in unique(samples_per_country_month$country_recoded)) {
    samples_per_month <- samples_per_country_month %>% 
      filter(country_recoded == country_c)
    for (y_m in samples_per_month$year_month) {
      available_samples <- metadata %>% 
        filter(country_recoded == country_c, year_month == y_m)
      n_samples <- samples_per_month[samples_per_month$year_month == y_m, "n_seqs_actual"]
      selected_samples_i <- available_samples %>%
        sample_n(size = n_samples, replace = F)
      if (is_first) {
        context_seq_data <- selected_samples_i
        is_first <- F
      } else {
        context_seq_data <- rbind(context_seq_data, selected_samples_i) 
      }
    }
  }
  
  context_seq_summary <- context_seq_data %>% 
    group_by(country_recoded, year_month) %>% 
    summarize(n_selected_samples = n())
  context_seq_summary <- merge(x = context_seq_summary, y = samples_per_country_month)
  
  write.table(
    x = context_seq_data,
    file = paste(OUTDIR, paste(PREFIX, "context_metadata.txt", sep = "_"), sep = "/"),
    row.names = F, col.names = T, quote = F, sep = "\t")
  
  selected_seq_strains <- unique(
    c(priority_seq_data$strain, context_seq_data$strain))
  selected_seq_data <- metadata[metadata$strain %in% selected_seq_strains, ]
  
  # Add outgroup seqs
  for (outgroup_strain in OUTGROUP_IDS) {
    if (!(outgroup_strain %in% selected_seq_data$strain)) {
      print(paste("Adding outroup seq", outgroup_strain))
      if (!(outgroup_strain %in% metadata$strain)) {
        stop(paste("Outgroup seq not in metadata."))
      }
      selected_seq_data_i <- metadata %>% filter(strain == outgroup_strain)
      selected_seq_data <- rbind(selected_seq_data, selected_seq_data_i)
    }
  }
  
  # Add swiss sequences
  selected_seq_data_swiss <- metadata %>% filter(country_recoded == "Switzerland")
  selected_seq_data <- rbind(selected_seq_data, selected_seq_data_swiss)
  
  # Change strain name to include GISAID_EPI_ISL and collection date
  selected_seq_data$old_strain <- selected_seq_data$strain
  selected_seq_data$strain <- paste(
    gsub(
      x = gsub(
        x = selected_seq_data$old_strain, pattern = "\\/", replacement = "_"),
      pattern = " ", replacement = "_"),
    selected_seq_data$gisaid_epi_isl,
    selected_seq_data$date, sep = "|")
  
  # Fetch selected sequences from the master alignment, write to file
  headers_found <- fetch_seqs_from_gisaid_alignment(
    gisaid_fp = ALIGNMENT,
    seq_headers = paste(">", selected_seq_data$old_strain, sep = ""),
    new_headers = paste(">", selected_seq_data$strain, sep = ""),
    seq_fp = seq_fp)

  # Check that all selected sequences were found
  not_found_filter <- !(
    paste(">", selected_seq_data$old_strain, sep = "") %in% headers_found)
  if (any(not_found_filter)) {
    print(selected_seq_data[not_found_filter, "old_strain"])
    warning(paste("These headers not found in", ALIGNMENT))
  }
  
  write.table(
    x = selected_seq_data,
    file = paste(OUTDIR, paste(PREFIX, "tree_metadata.txt", sep = "_"), sep = "/"),
    quote = F, row.names = F, col.names = T, sep = "\t")
}
