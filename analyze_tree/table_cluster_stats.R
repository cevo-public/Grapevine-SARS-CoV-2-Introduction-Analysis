require(dplyr)
require(lubridate)

# CLUSTERS <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-08-19_ch_cluster_analysis/clusters/rep_3_M_3_P_1_clusters.txt"
# METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-08-19_ch_cluster_analysis/data/n_context_1000_n_prority_1000/rep_3/tree_metadata.csv"

parser <- argparse::ArgumentParser()
parser$add_argument("--clusters", type="character", help="")
parser$add_argument("--metadata", type = "character", help = "")
parser$add_argument("--outdir", type = "character", help  = "")
parser$add_argument("--prefix", type = "character", help  = "")
parser$add_argument("--plotutilityfunctions", type = "character")
parser$add_argument("--nswissseqs", type = "integer")
args <- parser$parse_args()

CLUSTERS <- args$clusters
METADATA <- args$metadata
OUTDIR <- args$outdir
PREFIX <- args$prefix
UTILITY_FUNCTIONS <- args$plotutilityfunctions
N_SWISS_SEQS <- args$nswissseqs

source(UTILITY_FUNCTIONS)
system(command = paste("mkdir -p", OUTDIR))

# Load data
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
cluster_data <- read.delim(file = CLUSTERS)

cluster_data <- cluster_data %>%
  mutate(cluster_idx = 1:n())

cluster_data_by_tip <-get_cluster_data_by_tip(
  cluster_data = cluster_data, metadata = metadata)

cluster_summary <- cluster_data_by_tip %>%
  group_by(cluster_idx) %>%
  summarize(first_sample = min(`date`),
            last_sample = max(`date`))

# Add tmrca information to cluster summary
cluster_summary <- merge(
  x = cluster_summary, y = cluster_data, 
  by = c("cluster_idx"), all = T)

# Summarize how many sequences cluster and avg cluster size, sd
cluster_stats_summary <- cluster_summary %>%
  summarize(
    biggest_cluster = max(size),
    longest_lived_cluster = max(week(last_sample) - week(first_sample)),
    n_clusters = n(),
    mean_cluster_size = mean(size),
    sd_cluster_size = sd(size),
    n_singletons = sum(size == 1),
    n_transmission_chain = sum(size > 1),
    n_singletons = sum(size == 1),
    n_tips_excluded = N_SWISS_SEQS - sum(size),
    n_seqs_considered = nrow(metadata))

table_to_print <- cluster_stats_summary %>%
  select(biggest_cluster, longest_lived_cluster, n_transmission_chain, n_singletons, n_seqs_considered)

write.table(
  x = table_to_print,
  file = paste(OUTDIR, "/", PREFIX, "_cluster_stats.txt", sep = ""),
  quote = F, row.names = F, col.names = T, sep = "\t")


