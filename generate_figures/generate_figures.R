source("database/R/utility.R")
source("utility_functions.R")
source("generate_figures/functions.R")
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(ggtree)))
suppressMessages(suppressWarnings(require(argparse)))

# max_date <- "2020-12-01"
# workdir <- "/Users/nadeaus/Repos/grapevine/workdir"
# max_sampling_frac <- 0.05
# dates_to_highlight <- NULL
# min_chain_size <- 1
# focal_country <- "CHE"

parser <- argparse::ArgumentParser()
parser$add_argument("--maxdate", type="character")
parser$add_argument("--workdir", type="character")
parser$add_argument("--minchainsizeforplot", type="integer", default = 1)
parser$add_argument("--maxsamplingfrac", type="double", default = 1)
parser$add_argument("--datestohighlight", type="character", default = NULL)
parser$add_argument("--focalcountry", type="character")

args <- parser$parse_args()

max_date <- args$maxdate
workdir <- args$workdir
min_chain_size <- args$minchainsizeforplot
max_sampling_frac <- args$maxsamplingfrac
dates_to_highlight <- args$datestohighlight
focal_country <- args$focalcountry
if (!is.null(dates_to_highlight)) {
  dates_to_highlight <- strsplit(dates_to_highlight, split = " ")[[1]]
}

db_connection = open_database_connection()
outdir <- paste(workdir, "output", sep = "/")

print("Defining variables to be shared across figures.")
country_colors <- get_country_colors(db_connection)

print("Plotting barchart of sampling intensity through time.")
sampling_intensity_data <- plot_sampling_intensity(
  db_connection = db_connection,
  workdir = workdir,
  outdir = outdir,
  max_date = max_date,
  max_sampling_frac = max_sampling_frac,
  dates_to_highlight = dates_to_highlight,
  focal_country = focal_country
)

est_imports_file <- paste(workdir, "tmp/alignments/estimated_travel_cases.csv", sep = "/")
if (file.exists(est_imports_file)) {
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
} else {
  print(paste(
    "Not plotting transmission chain origins through time because file",
    est_imports_file,"not found."))
}

print("Plotting transmission chains through time.")
plot_chains(
  workdir = workdir,
  outdir = outdir,
  country_colors = country_colors,
  min_chain_size = min_chain_size,
  plot_height_in = 15
)

print("Plotting number of introductions, chain extinctions through time.")
p <- plot_introductions_and_extinctions(
  workdir = workdir,
  outdir = outdir
)



