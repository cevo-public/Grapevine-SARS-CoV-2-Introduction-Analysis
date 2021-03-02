install.packages(c(
  "argparse",
  "lubridate",
  "scatterpie",
  "tidyr",
  "tidytree",
  "dplyr",
  "config",
  "DBI",
  "countrycode",
  "RColorBrewer"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("treeio")
BiocManager::install("ggtree")