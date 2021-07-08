source("database/R/utility.R")
source("generate_figures/functions.R")
source("utility_functions.R")
require(dplyr)

# workdir <- "/Users/nadeaus/Repos/grapevine/workdir"
# max_date <- "2020-11-30"
# focal_country <- "NZL"

parser <- argparse::ArgumentParser()
parser$add_argument("--workdir", type = "character")
parser$add_argument("--maxdate", type = "character")
parser$add_argument("--focalcountry", type = "character")

args <- parser$parse_args()

workdir <- args$workdir
max_date <- as.Date(args$maxdate)
focal_country <- as.Date(args$focalcountry)

db_connection = open_database_connection()
outdir <- paste(workdir, "output/transmission_chain_alignments", sep = "/")
system(command = paste("mkdir -p", outdir))

print("Loading sample metadata and inferred transmission chain data.")
sample_metadata <- load_sample_metadata(workdir = workdir)
chains_max <- load_chain_asr_data(l = F, workdir = workdir) %>%
  mutate(chain_idx = 1:n(),
         chains_assumption = "max")
chains_min <- load_chain_asr_data(l = T, workdir = workdir) %>%
  mutate(chain_idx = 1:n(),
         chains_assumption = "min")

print(paste("Filtering samples for BDSKY analysis to dates <=", max_date))
samples_to_write_out <- rbind(
  pivot_chains_to_samples(chains = chains_max, metadata = sample_metadata),
  pivot_chains_to_samples(chains = chains_min, metadata = sample_metadata)) %>%
  tidyr::unite(col = "header", strain, sample, chain_idx, sep = "|", remove = F) %>%  # headers form: <strain with "\" & " " --> "_">|<gisaid_epi_isl>|<yyyy-mm-dd date>|<cluster_idx>
  filter(date <= as.Date(max_date), focal_sequence)

if (focal_country == "CHE") {
  print("Filtering sample metadata to Viollier only, adding transmission chain information.")
  samples_to_write_out <- samples_to_write_out %>%
    filter(originating_lab == "Viollier AG",
           submitting_lab == "Department of Biosystems Science and Engineering, ETH Zürich")  # because some samples sequenced by University Hospital Basel, Clinical Bacteriology also come from Viollier
}
  
max_chain_header_mapping <- unlist(
  samples_to_write_out %>%
    filter(chains_assumption == "max") %>%
    select(header))
names(max_chain_header_mapping) <- unlist(
  samples_to_write_out %>%
    filter(chains_assumption == "max") %>%
    select(gisaid_epi_isl))
min_chain_header_mapping <- unlist(
  samples_to_write_out %>%
    filter(chains_assumption == "min") %>%
    select(header))
names(min_chain_header_mapping) <- unlist(
  samples_to_write_out %>%
    filter(chains_assumption == "min") %>%
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

print("Writing out summary information for bdsky setup.")
alignment_summary <- samples_to_write_out %>%
  group_by(chains_assumption) %>%
  summarise(n_seqs_in_alignment = n(),
            n_unique_chains = length(unique(chain_idx)))
chain_summary <- rbind(chains_max, chains_min) %>% 
  group_by(chains_assumption) %>%
  summarise(n_singleton_chains = sum(size == 1),
            n_seqs_in_biggest_chain = max(size))
summary <- merge(alignment_summary, chain_summary)

write.csv(x = summary,
          file = paste(outdir, "chains_summary.csv", sep = "/"),
          row.names = F)
