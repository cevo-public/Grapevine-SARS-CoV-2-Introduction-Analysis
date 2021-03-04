source("database/R/utility.R")
source("utility_functions.R")
source("generate_figures/functions.R")
require(dplyr)
require(ggplot2)
require(ggtree)
require(argparse)

# min_date <- "2020-01-01"
# max_date <- "2020-12-31"
# min_length <- 27000
# workdir <- "/Users/nadeaus/NonRepoProjects/cov-swiss-phylogenetics/grapevine/jan-dec_no-max-sampling_-5_context-sf"

parser <- argparse::ArgumentParser()
parser$add_argument("--mindate", type="character")
parser$add_argument("--maxdate", type="character")
parser$add_argument("--minlength", type="integer")
parser$add_argument("--workdir", type="character")

args <- parser$parse_args()

min_date <- args$mindate
max_date <- args$maxdate
min_length <- args$minlength
workdir <- args$workdir

db_connection = open_database_connection()
outdir <- paste(workdir, "output", sep = "/")
system(command = paste("mkdir -p", outdir))

qcd_gisaid_query <- dplyr::tbl(db_connection, "gisaid_sequence") %>%
  filter(
    date <= !! max_date,
    date >= !! min_date, 
    length >= min_length,
    host == 'Human', 
    nextclade_qc_snp_clusters_status == 'good', 
    nextclade_qc_private_mutations_status == 'good', 
    nextclade_qc_overall_status != 'bad')

# Define variables to be shared across figures
country_colors <- get_country_colors(db_connection)

# # Plot barchart of sampling intensity through time
# sampling_intensity_data <- plot_sampling_intensity(
#   db_connection = db_connection,
#   qcd_gisaid_query = qcd_gisaid_query,
#   outdir = outdir
# )
#
# # Plot barchart of transmission chain origins through time
# for (s in c(T, F)) {
#   plot_chain_sources(
#     s = s,
#     workdir = workdir,
#     outdir = outdir,
#     country_colors = country_colors
#   )
# }

# Plot transmission chains through time
plot_chains(
  workdir = workdir,
  outdir = outdir,
  country_colors = country_colors,
  min_chain_size = 3,
  plot_height_in = 15
)






