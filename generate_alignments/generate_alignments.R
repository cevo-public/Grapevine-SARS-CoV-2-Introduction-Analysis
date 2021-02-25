source("database/R/utility.R")
source("utility_functions.R")
source("generate_alignments/functions.R")
require(dplyr)

# min_date <- "2020-08-01"
# max_date <- "2020-08-31"
# min_length <- 27000
# travel_context_scale_factor <- 5
# similarity_context_scale_factor <- 1
# outdir <- "/Users/nadeaus/Repos/grapevine/dont_commit/aug_test/alignments"
# python_path <- "/Users/nadeaus/Repos/database/python/venv/bin/python3"
# reference <- "/Users/nadeaus/Repos/database/python/ncov/defaults/reference_seq.fasta"

parser <- argparse::ArgumentParser()
parser$add_argument("--mindate", type="character")
parser$add_argument("--maxdate", type="character")
parser$add_argument("--minlength", type="integer")
parser$add_argument("--travelcontextscalefactor", type="integer", help="Multiplicative factor, how many times the # swiss sequences should we select for the travel context set?")
parser$add_argument("--similaritycontextscalefactor", type="integer", help="Multiplicative factor, how many times the # swiss sequences should we select for the genetic similarity context set?")
parser$add_argument("--outdir", type="character")
parser$add_argument("--pythonpath", type="character", help="e.g. /Users/nadeaus/Repos/database/python/venv/bin/python3")
parser$add_argument("--reference", type="character", help="e.g. /Users/nadeaus/Repos/database/python/ncov/defaults/reference_seq.fasta")

args <- parser$parse_args()

min_date <- args$mindate
max_date <- args$maxdate
min_length <- args$minlength
travel_context_scale_factor <- args$travelcontextscalefactor
similarity_context_scale_factor <- args$similaritycontextscalefactor
outdir <- args$outdir
python_path <- args$pythonpath
reference <- args$reference

# Hardcoded parameters
outgroup_gisaid_epi_isls = c("EPI_ISL_406798", "EPI_ISL_402125")  # The nextstrain global tree is rooted between these two sequences (Wuhan/WH01/2019 & Wuhan/Hu-1/2019), which you can see by filtering the tree to Chinese sequences (to make it reasonably small), downloading the newick tree, and plotting it.

db_connection = open_database_connection()
system(command = paste("mkdir -p", outdir))

# QC GISAID sequences
qcd_gisaid_query <- dplyr::tbl(db_connection, "gisaid_sequence") %>%
  filter(
    date <= !! max_date,
    date >= !! min_date, 
    length >= min_length,
    host == 'Human', 
    nextclade_qc_snp_clusters_status == 'good', 
    nextclade_qc_private_mutations_status == 'good', 
    nextclade_qc_overall_status != 'bad')

# Get pangolin lineages to write out alignments for
lineages <- get_pangolin_lineages(
  db_connection = db_connection,
  outdir = outdir,
  qcd_gisaid_query = qcd_gisaid_query
)

# Estimate travel cases by source country and month
travel_cases <- get_travel_cases(
  db_connection = db_connection,
  min_date = min_date,
  max_date = max_date,
  outdir = outdir
)

# Get travel context set based on travel cases
travel_strains <- get_travel_strains(
  n_strains = sum(lineages$is_swiss_TRUE) * travel_context_scale_factor,
  travel_cases = travel_cases,
  db_connection = db_connection,
  outdir = outdir,
  qcd_gisaid_query = qcd_gisaid_query
)

# Get similarity context set based on genetic proximity
similarity_strains <- get_similarity_strains(
  similarity_context_scale_factor = similarity_context_scale_factor,
  lineages = lineages,
  db_connection = db_connection,
  qcd_gisaid_query = qcd_gisaid_query,
  python_path = python_path,
  reference = reference
)

# Write out alignments, one per pangolin lineage
alignments <- write_out_alignments(
  lineages = lineages,
  travel_strains = travel_strains,
  similarity_strains = similarity_strains,
  outgroup_gisaid_epi_isls = outgroup_gisaid_epi_isls,
  outdir = outdir,
  qcd_gisaid_query = qcd_gisaid_query,
  db_connection = db_connection
)

write.csv(
  x = alignments,
  file = paste(outdir, "alignment_sizes.csv", sep = "/"),
  row.names = F)





