source("database/R/utility.R")
source("utility_functions.R")
source("figures/functions.R")
require(dplyr)
require(ggplot2)
require(ggtree)

min_date <- "2020-01-01"
max_date <- "2021-03-01"
min_length <- 27000

db_connection = open_database_connection()
workdir <- "/Users/nadeaus/Repos/grapevine/workdir"
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

# Plot barchart of sampling intensity through time
sampling_intensity_data <- plot_sampling_intensity(
  db_connection = db_connection,
  qcd_gisaid_query = qcd_gisaid_query,
  outdir = outdir
)

# Plot barchart of transmission chain origins through time
for (s in c(T, F)) {
  plot_chain_sources(
    s = s, 
    workdir = workdir,
    outdir = outdir, 
    country_colors = country_colors
  )
}






