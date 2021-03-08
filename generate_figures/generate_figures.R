source("database/R/utility.R")
source("utility_functions.R")
source("generate_figures/functions.R")
require(dplyr)
require(ggplot2)
require(ggtree)
require(argparse)

# max_date <- "2020-12-31"
# workdir <- "/Users/nadeaus/NonRepoProjects/cov-swiss-phylogenetics/grapevine/jan-dec_-005_max-sampling_-5_context-sf"

parser <- argparse::ArgumentParser()
parser$add_argument("--maxdate", type="character")
parser$add_argument("--workdir", type="character")

args <- parser$parse_args()

min_date <- args$mindate
max_date <- args$maxdate
min_length <- args$minlength
workdir <- args$workdir

db_connection = open_database_connection()
outdir <- paste(workdir, "output", sep = "/")
system(command = paste("mkdir -p", outdir))

# Define variables to be shared across figures
country_colors <- get_country_colors(db_connection)

# Plot barchart of sampling intensity through time
sampling_intensity_data <- plot_sampling_intensity(
  db_connection = db_connection,
  workdir = workdir,
  outdir = outdir,
  max_date = max_date
)

# Plot estimated infectious arrivals through time
plot_source_prior(
  workdir = workdir,
  outdir = outdir,
  country_colors = country_colors
)

# Plot transmission chain origins through time
for (s in c(T, F)) {
  plot_chain_sources(
    s = s,
    workdir = workdir,
    outdir = outdir,
    country_colors = country_colors
  )
}

# Plot transmission chains through time
plot_chains(
  workdir = workdir,
  outdir = outdir,
  country_colors = country_colors,
  min_chain_size = 3,
  plot_height_in = 15
)

# Plot number of introductions, chain deaths through time
plot_introductions_and_extinctions(
  workdir = workdir,
  outdir = outdir
)

# Plot introductions/time during on-qurantine list vs. off-quarantine list periods
# TODO
plot_travel_quarantine_effect(
  workdir = workdir,
  outdir = outdir,
  db_connection = db_connection
)




