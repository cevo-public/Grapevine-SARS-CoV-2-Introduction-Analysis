# Create date to date_trun week lookup table for Tim.
# Will allow us to have the same Re breakpoints as sampling breakpoints.

source("database/R/utility.R")
source("utility_functions.R")
suppressMessages(suppressWarnings(require(dplyr)))

# outdir <- "/Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/jan-dec_-01_max_sampling_1_travel_1_sim_context-sf_111_travel-wt/output/transmission_chain_alignments"

parser <- argparse::ArgumentParser()
parser$add_argument("--outdir", type = "character")
parser$add_argument("--maxdate", type = "character")
parser$add_argument("--mindate", type = "character", default = "2020-01-01")

args <- parser$parse_args()

outdir <- args$outdir
max_date <- as.Date(args$maxdate)
min_date <- as.Date(args$mindate)

# Connect to database
db_connection = open_database_connection()

# Generate date data
dates <- data.frame(date = seq.Date(from = min_date, to = max_date, by = "day"))

# Create temporary table for date mapping
tbl_name <- "date_to_week"
DBI::dbCreateTable(
  conn = db_connection,
  name = tbl_name,
  fields = dates
)

DBI::dbWriteTable(
  conn = db_connection,
  name = tbl_name, 
  value = dates,
  overwrite = T
)

# Generate mapping
dates_with_week <- dplyr::tbl(db_connection, tbl_name) %>%
  mutate(week = date_trunc('week', date)) %>%
  collect() %>%
  mutate(week_orig = week,
         week = format(week_orig, "%Y-%m-%d"))

# Write out mapping
write.csv(
  x = dates_with_week %>% select(date, week),
  file = paste(outdir, "date_to_week.csv", sep = "/"),
  row.names = F
)

# Remove table
DBI::dbRemoveTable(
  conn = db_connection,
  name = tbl_name
)
