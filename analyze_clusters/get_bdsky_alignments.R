source("database/R/utility.R")
source("generate_figures/functions.R")
source("utility_functions.R")
require(dplyr)

# workdir <- "/Users/nadeaus/NonRepoProjects/cov-swiss-phylogenetics/grapevine/jan-dec_-005_max-sampling_-5_context-sf"

parser <- argparse::ArgumentParser()
parser$add_argument("--workdir", type = "character")
parser$add_argument("--maxdate", type = "character")

args <- parser$parse_args()

workdir <- args$workdir
max_date <- as.Date(args$maxdate)

db_connection = open_database_connection()
outdir <- paste(workdir, "output/transmission_chain_alignments", sep = "/")
system(command = paste("mkdir -p", outdir))

print("Loading sample metadata and inferred transmission chain data.")
sample_metadata <- load_sample_metadata(workdir = workdir)
chains_max <- load_chain_asr_data(s = F, workdir = workdir) %>%
  mutate(chain_idx = 1:n(),
         chains_assumption = "max")
chains_min <- load_chain_asr_data(s = T, workdir = workdir) %>%
  mutate(chain_idx = 1:n(),
         chains_assumption = "min")

print("Filtering sample metadata to Viollier only, adding transmission chain information.")
viollier_samples <- rbind(
  pivot_chains_to_samples(chains = chains_max, metadata = sample_metadata),
  pivot_chains_to_samples(chains = chains_min, metadata = sample_metadata)) %>%
  filter(originating_lab == "Viollier AG", 
         submitting_lab == "Department of Biosystems Science and Engineering, ETH Zürich",  # because some samples sequenced by University Hospital Basel, Clinical Bacteriology also come from Viollier
         focal_sequence) %>%  
  tidyr::unite(col = "header", strain, sample, chain_idx, sep = "|", remove = F) %>%  # headers form: <strain with "\" & " " --> "_">|<gisaid_epi_isl>|<yyyy-mm-dd date>|<cluster_idx>
  mutate(after_may_1 = date >= as.Date("2020-05-01")) %>%
  filter(date <= max_date)
  
max_chain_header_mapping <- unlist(
  viollier_samples %>%
    filter(chains_assumption == "max") %>%
    select(header))
names(max_chain_header_mapping) <- unlist(
  viollier_samples %>%
    filter(chains_assumption == "max") %>%
    select(gisaid_epi_isl))
min_chain_header_mapping <- unlist(
  viollier_samples %>%
    filter(chains_assumption == "min") %>%
    select(header))
names(min_chain_header_mapping) <- unlist(
  viollier_samples %>%
    filter(chains_assumption == "min") %>%
    select(gisaid_epi_isl))

max_chain_after_may_1_header_mapping <- unlist(
  viollier_samples %>%
    filter(chains_assumption == "max", after_may_1) %>%
    select(header))
names(max_chain_after_may_1_header_mapping) <- unlist(
  viollier_samples %>%
    filter(chains_assumption == "max", after_may_1) %>%
    select(gisaid_epi_isl))
min_chain_after_may_1_header_mapping <- unlist(
  viollier_samples %>%
    filter(chains_assumption == "min", after_may_1) %>%
    select(header))
names(min_chain_after_may_1_header_mapping) <- unlist(
  viollier_samples %>%
    filter(chains_assumption == "min", after_may_1) %>%
    select(gisaid_epi_isl))

print("Writing out alignments.")
export_seqs_as_fasta(
  db_connection = db_connection,
  sample_names = names(max_chain_header_mapping),
  seq_outfile = paste(outdir, "max_chains.fasta", sep = "/"),
  table = "gisaid_sequence",
  sample_name_col = "gisaid_epi_isl",
  header_mapping = max_chain_header_mapping,
  seq_col = "aligned_seq")
export_seqs_as_fasta(
  db_connection = db_connection,
  sample_names = names(min_chain_header_mapping),
  seq_outfile = paste(outdir, "min_chains.fasta", sep = "/"),
  table = "gisaid_sequence",
  sample_name_col = "gisaid_epi_isl",
  header_mapping = min_chain_header_mapping,
  seq_col = "aligned_seq")

# Repeat analysis with only samples beginning May 1
export_seqs_as_fasta(
  db_connection = db_connection,
  sample_names = names(max_chain_after_may_1_header_mapping),
  seq_outfile = paste(outdir, "max_chains_after_may_1.fasta", sep = "/"),
  table = "gisaid_sequence",
  sample_name_col = "gisaid_epi_isl",
  header_mapping = max_chain_after_may_1_header_mapping,
  seq_col = "aligned_seq")
export_seqs_as_fasta(
  db_connection = db_connection,
  sample_names = names(min_chain_after_may_1_header_mapping),
  seq_outfile = paste(outdir, "min_chains_after_may_1.fasta", sep = "/"),
  table = "gisaid_sequence",
  sample_name_col = "gisaid_epi_isl",
  header_mapping = min_chain_after_may_1_header_mapping,
  seq_col = "aligned_seq")

print("Writing out summary information for bdsky setup.")
alignment_summary <- viollier_samples %>% 
  group_by(chains_assumption) %>%
  summarise(n_seqs_in_alignment = n(),
            n_unique_chains = length(unique(chain_idx)))
chain_summary <- rbind(chains_max, chains_min) %>% 
  group_by(chains_assumption) %>%
  summarise(n_singleton_chains = sum(size == 1),
            n_seqs_in_biggest_chain = max(size))
summary <- merge(alignment_summary, chain_summary)  %>%
  mutate(timespan = "after jan 1")

alignment_after_may_1_summary <- viollier_samples %>% 
  filter(after_may_1) %>%
  group_by(chains_assumption) %>%
  summarise(n_seqs_in_alignment = n(),
            n_unique_chains = length(unique(chain_idx)))
chain_after_may_1_summary <- viollier_samples %>%
  filter(after_may_1) %>%
  group_by(chain_idx, chains_assumption) %>%
  summarise(size = n()) %>%
  ungroup() %>%
  group_by(chains_assumption) %>%
  summarise(n_singleton_chains = sum(size == 1),
            n_seqs_in_biggest_chain = max(size))

summary_after_may_1 <- merge(alignment_after_may_1_summary, chain_after_may_1_summary) %>%
  mutate(timespan = "after may 1")

summary_to_print <- rbind(x = summary, y = summary_after_may_1)

write.csv(x = summary_to_print,
          file = paste(outdir, "chains_summary.csv", sep = "/"),
          row.names = F)
