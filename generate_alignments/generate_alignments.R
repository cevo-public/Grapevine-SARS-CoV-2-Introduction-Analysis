source("database/R/utility.R")
source("utility_functions.R")
source("generate_alignments/functions.R")
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(yaml)))
suppressMessages(suppressWarnings(require(stringr)))

# config <- "/Users/nadeaus/Repos/grapevine/example_workdir/input/grapevine_config.yaml"
# outdir <- "~/Downloads"
# python_path <- "/Users/nadeaus/Repos/database/python/venv/bin/python3"
# reference <- "/Users/nadeaus/Repos/database/python/ncov/defaults/reference_seq.fasta"

parser <- argparse::ArgumentParser()
parser$add_argument("--config", type="character", help = "path to grapevine_config.yaml file.")
parser$add_argument("--outdir", type="character")
parser$add_argument("--pythonpath", type="character", help="Path to python3 with required packages installed.")
parser$add_argument("--reference", type="character", help="Reference sequence.")

args <- parser$parse_args()

config <- args$config
outdir <- args$outdir
python_path <- args$pythonpath
reference <- args$reference

config_values <- yaml::read_yaml(file = config)
focal_country <- config_values$focal_country
max_date <- config_values$max_date
min_date <- config_values$min_date
max_missing <- config_values$max_missing
min_length <- config_values$min_length
max_sampling_frac <- config_values$max_sampling_fraction
subsample_by_canton <- config_values$subsample_by_canton
smooth_conf_cases <- config_values$smooth_conf_cases
travel_context_scale_factor <- config_values$travel_context_scale_factor
similarity_context_scale_factor <- config_values$similarity_context_scale_factor
travel_data_weights <- config_values$travel_data_weights
which_trees <- config_values$which_trees
n_trees <- config_values$n_trees
outgroup_gisaid_epi_isls <- strsplit(config_values$outgroup_gisaid_epi_isls, split = " ")[[1]]
unique_context_only <- config_values$unique_context_only
favor_exposures <- config_values$favor_exposures
mask_from_start <- config_values$mask_from_start
mask_from_end <- config_values$mask_from_end

db_connection = open_database_connection()
system(command = paste("mkdir -p", outdir))

# QC GISAID sequences
qcd_gisaid_query <- dplyr::tbl(db_connection, "gisaid_api_sequence") %>%
  filter(
    date <= !! max_date,
    date >= !! min_date,
    !is.na(date),
    nextclade_total_missing <= max_missing,
    nextclade_alignment_end - nextclade_alignment_start + 1 >= min_length,
    host == 'Human',
    nextclade_qc_snp_clusters_status == 'good',
    nextclade_qc_private_mutations_status == 'good',
    nextclade_qc_overall_status != 'bad')

if (which_trees != '\\.*') {
  qcd_gisaid_query <- qcd_gisaid_query %>%
  filter(grepl(x = pangolin_lineage, pattern = which_trees))
  unique_lineages <- qcd_gisaid_query %>%
    distinct(pangolin_lineage) %>%
    collect()
  cat(paste0(
    "Specified which_trees = '", which_trees, "' ",
    "so not including nextclade qc filters and", 
    " only including sequences from these lineages:\n", 
    paste0(unique_lineages$pangolin_lineage, collapse = "\n"),
    "\n"))
}

# Get pangolin lineages to write out alignments for
lineages <- get_pangolin_lineages(
  db_connection = db_connection,
  outdir = outdir,
  qcd_gisaid_query = qcd_gisaid_query,
  focal_country = focal_country
)

# Select focal sequences
qcd_gisaid_query <- select_sequences(
  qcd_gisaid_query = qcd_gisaid_query,
  max_sampling_frac = max_sampling_frac,
  focal_country = focal_country,
  db_connection = db_connection,
  outdir = outdir,
  subsample_by_canton = subsample_by_canton,
  smooth_conf_cases = smooth_conf_cases,
  max_date = max_date,
  min_date = min_date
)

if (n_trees > 0) {
  n_lineages <- min(n_trees, nrow(lineages))
  warning(paste("Only generating", n_lineages, "trees even though there are", nrow(lineages), "total lineages."))
  lineages <- lineages %>% arrange(is_swiss_TRUE)
  lineages <- lineages[1:n_lineages, ]
}

# Select location context sequences
if (travel_context_scale_factor > 0 & focal_country == "CHE") {
  # Estimate travel cases by source country and month
  travel_cases <- get_travel_cases(
    db_connection = db_connection,
    min_date = min_date,
    max_date = max_date,
    travel_data_weights = travel_data_weights,
    outdir = outdir
  )

  # Get travel context set based on travel cases
  n_strains <- nrow(qcd_gisaid_query %>%
    filter(iso_country == "CHE") %>%
    collect())
  travel_strains <- get_travel_strains(
    n_strains = ceiling(n_strains * travel_context_scale_factor),
    travel_cases = travel_cases,
    db_connection = db_connection,
    outdir = outdir,
    qcd_gisaid_query = qcd_gisaid_query
  )
} else {
  # Return dummy travel strains
  travel_strains <- data.frame(
      gisaid_epi_isl = NA,
      iso_country = NA,
      date = NA,
      date_str = NA,
      pangolin_lineage = NA)
}

# Select genetic proximity context sequences
similarity_strains <- get_similarity_strains(
  similarity_context_scale_factor = similarity_context_scale_factor,
  unique_context_only = unique_context_only,
  lineages = lineages,
  db_connection = db_connection,
  qcd_gisaid_query = qcd_gisaid_query,
  python_path = python_path,
  reference = reference,
  focal_country = focal_country
)

# Write out alignments, one per pangolin lineage
alignments <- write_out_alignments(
  lineages = lineages,
  travel_strains = travel_strains,
  similarity_strains = similarity_strains,
  outgroup_gisaid_epi_isls = outgroup_gisaid_epi_isls,
  outdir = outdir,
  qcd_gisaid_query = qcd_gisaid_query,
  db_connection = db_connection,
  mask_from_start = mask_from_start,
  mask_from_end = mask_from_end,
  focal_country = focal_country
)

write.csv(
  x = alignments,
  file = paste(outdir, "alignment_sizes.csv", sep = "/"),
  row.names = F)





