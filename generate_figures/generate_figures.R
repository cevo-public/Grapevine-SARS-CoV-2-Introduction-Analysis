source("database/R/utility.R")
source("utility_functions.R")
source("figures/functions.R")
require(dplyr)
require(ggplot2)
require(ggtree)

db_connection = open_database_connection()
workdir <- "/Users/nadeaus/Repos/grapevine/dont_commit/test_travel_scale_2"
outdir <- paste(workdir, "output", sep = "/")

country_colors <- get_country_colors(db_connection)

for (s in c(T, F)) {
  plot_chain_sources(
    s = s, 
    workdir = workdir,
    outdir = outdir, 
    country_colors = country_colors
  )
}


