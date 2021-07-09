# test pangolin lineage rollup

source("generate_alignments/functions.R")
source("utility_functions.R")
db_connection <- open_database_connection("server")

raw_data <- read.csv("database/data/pangolin_lineage_table.csv")

clean_data <- raw_data %>%
  filter(grepl(x = Description, pattern = "(A|a)lias"),
         !(grepl(x = Description, pattern = "(W|w)ithdrawn"))) %>%
  mutate(
    alias = stringr::str_extract(
      string = Description, 
      pattern = " [:alpha:](\\.[:digit:]+)*( |,|$)"),
    alias = gsub(x = alias, pattern = "( |,)", replacement = ""))

alias_table <- dplyr::tbl(db_connection, "pangolin_lineage_alias") %>% collect()

results <- clean_data %>%
  select(Lineage, alias) %>%
  mutate(
    guess = unlist(lapply(
      FUN = expand_pangolin_lineage_aliases,
      X = Lineage,
      alias_table = alias_table)),
    parent_guess = unlist(lapply(
      FUN = get_parent_lineage,
      X = guess)),
    grandparent_guess = unlist(lapply(
      FUN = get_parent_lineage,
      X = parent_guess)),
    greatgrandparent_guess = unlist(lapply(
      FUN = get_parent_lineage,
      X = grandparent_guess)),
    greatgreatgrandparent_guess = unlist(lapply(
      FUN = get_parent_lineage,
      X = greatgrandparent_guess)))

nrow(results %>% filter(alias != guess)) == 0

