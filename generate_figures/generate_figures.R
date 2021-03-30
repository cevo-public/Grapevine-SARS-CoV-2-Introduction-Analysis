source("database/R/utility.R")
source("utility_functions.R")
source("generate_figures/functions.R")
require(dplyr)
require(ggplot2)
require(ggtree)
require(argparse)

max_date <- "2020-12-31"
workdir <- "/Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/jan-dec_-01_max_sampling_1_travel_-5_sim_context-sf_1_exp-wt"
min_chain_size <- 1

# parser <- argparse::ArgumentParser()
# parser$add_argument("--maxdate", type="character")
# parser$add_argument("--workdir", type="character")
# parser$add_argument("--minchainsizeforplot", type="integer", default = 1)
# 
# args <- parser$parse_args()
# 
# max_date <- args$maxdate
# workdir <- args$workdir
# min_chain_size <- args$minchainsizeforplot

db_connection = open_database_connection()
outdir <- paste(workdir, "output", sep = "/")

print("Defining variables to be shared across figures.")
country_colors <- get_country_colors(db_connection)

print("Plotting barchart of sampling intensity through time.")
sampling_intensity_data <- plot_sampling_intensity(
  db_connection = db_connection,
  workdir = workdir,
  outdir = outdir,
  max_date = max_date
)

grapevine_results <- load_grapevine_results(
  workdir = workdir, 
  min_chain_size = 1,
  viollier_only = F)

print("Plotting transmission chain origins through time.")
plot_chain_origins(
  workdir = workdir,
  grapevine_results = grapevine_results,
  outdir = outdir,
  country_colors = country_colors)

print("Tabling prior vs. posterior transmission chain origins in 1st & 2nd wave.")
for (s in c(T, F)) {
  table_chain_origins(
    s = s,
    workdir = workdir,
    outdir = outdir
  )
}

print("Plotting transmission chains through time.")
plot_chains(
  workdir = workdir,
  outdir = outdir,
  country_colors = country_colors,
  min_chain_size = min_chain_size,
  plot_height_in = 15
)

print("Plotting number of introductions, chain extinctins through time.")
plot_introductions_and_extinctions(
  workdir = workdir,
  outdir = outdir
)



