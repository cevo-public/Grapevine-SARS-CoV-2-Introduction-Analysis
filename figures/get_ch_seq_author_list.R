# Get GISAID acknowledgments for Swiss sequnces

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-29_ch_cluster_analysis"
METADATA <- paste(WORKDIR, "/data/est_imports/metadata_all.txt", sep = "")
SWISS_SEQS <- paste(WORKDIR, "/data/qc_master_alignment/swiss_alignment_filtered2_masked_oneline.fasta", sep = "")

swiss_seqs <- ape::read.FASTA(file = SWISS_SEQS)
metadata <- read.delim(file = METADATA, stringsAsFactors = F, sep = "\t", quote = "")
ch_metadata <- metadata %>% filter(strain %in% names(swiss_seqs))

ch_author_info <- ch_metadata %>% 
  group_by(authors, paper_url) %>%
  summarise(n_seqs = n()) %>%
  arrange(desc(n_seqs))
View(ch_author_info)

# Most recent sequence from CH
View(ch_metadata %>% arrange(desc(date)) %>% head(5))
       