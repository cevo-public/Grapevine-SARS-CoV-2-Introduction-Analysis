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
  "RColorBrewer",
  "RPostgres",
  "gplots",
  "ggpubr"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

install.packages("devtools")
devtools::install_github("GuangchuangYu/treeio")