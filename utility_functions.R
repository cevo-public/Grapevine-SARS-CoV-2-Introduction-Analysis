# Utility functions for project, mostly having to do with manipulating data and 
# trees.

#' Overwrites the open_database_connection function: Connection data (including the password) must be provided via
#' environment variables.
# open_database_connection <- function (...) {
#   db_connection <- DBI::dbConnect(
#     RPostgres::Postgres(),
#     host = Sys.getenv("DB_HOST"),
#     port = Sys.getenv("DB_PORT"),
#     user = Sys.getenv("DB_USER"),
#     password = Sys.getenv("DB_PASSWORD"),
#     dbname = Sys.getenv("DB_DBNAME")
#   )
#   return(db_connection)
# }
open_database_connection <- function (
  db_instance = "server",
  config_file = "workdir/input/config.yml"
) {
  print(paste("config_file location:", config_file))
  print(paste("config file exists", file.exists(config_file)))
  connection_data <- config::get("database", file = config_file)[[db_instance]]
  db_connection <- DBI::dbConnect(
    RPostgres::Postgres(),
    host = connection_data$host,
    port = connection_data$port,
    user = connection_data$username,
    password = connection_data$dbpassword,
    dbname = connection_data$dbname
  )
  return(db_connection)
}



is_tip <- function(node, n_tips) {
  # Return T if the node is a tip based on the default node #ing and # tips
  return(as.numeric(node) <= n_tips)
}

#' Get node data for children of node
#' @param node integer node number
#' @param tree_data Tree data in tidytree treedata format.
#' @return rows of tree_data corresponding to the children of node.
get_child_node_data <- function(node, tree_data) {
  child_node_data <- tree_data[tree_data$parent == node, ]
  child_node_data <- child_node_data[child_node_data$parent != child_node_data$node, ]  # If node is root, do not return root as child
  return(child_node_data)
}

get_n_descendent_swiss_mrcas <- function(node, tree_data, chains, n_tips) {
  # Returns 1 if the node is a Swiss MRCA, otherwise returns the # Swiss MRCAs
  # descending from the node. Counts only the uppermost MRCA if there additional
  # Swiss MRCAs are nested under a Swiss MRCA higher in the tree.
  if (node %in% chains$ch_mrca) {
    return(1)
  } else if (node <= n_tips) {
    return(0)
  } else {
    child_node_data <- get_child_node_data(node = node, tree_data = tree_data)
  }
  n_descendent_swiss_mrcas <- sum(unlist(lapply(
    X = child_node_data$node,
    FUN = get_n_descendent_swiss_mrcas,
    tree_data = tree_data, 
    chains = chains, 
    n_tips = n_tips)))
  return(n_descendent_swiss_mrcas)
}

collapse_list <- function(list) {
  return(paste0(list, collapse = ", "))
}

get_n_context_children <- function(node, tree_data, context_node_set, n_tips) {
  child_node_data <- get_child_node_data(node, tree_data)
  child_has_context_seqs <- unlist(lapply(
    X = child_node_data$node, 
    FUN = get_node_has_context_seqs,
    tree_data = tree_data,
    context_node_set = context_node_set,
    n_tips = n_tips))
  return(sum(child_has_context_seqs))
}

recode_colnames <- function(col_name) {
  # Recode ASR locations (format "First.Second_loc_weight") to match country_recoded names
  if (grepl(x = col_name, pattern = "_loc_weight")) {
    col_name_new <- gsub(x = col_name, pattern = "_loc_weight", replacement = "")
    col_name_new <- gsub(x = col_name_new, pattern = "\\.", replacement = " ")
    return(col_name_new)
  }
  return(col_name)
}

get_node_has_context_seqs <- function(tree_data, node, context_node_set, n_tips) {
  # Return whether the node or any of its descendents are in the context node set
  is_tip <- node <= n_tips
  has_context <- node %in% context_node_set 
  while (!(is_tip | has_context)) {
    child_node_data <- get_child_node_data(node = node, tree_data = tree_data)
    has_context <- any(unlist(lapply(
      X = child_node_data$node, FUN = get_node_has_context_seqs, 
      tree_data = tree_data, context_node_set = context_node_set, n_tips = n_tips)))
  }
  return(has_context)
}

check_all_tips_in_metadata <- function(metadata_ids, tree_ids) {
  # Check all samples in tree present in metadata
  missing_in_metadata <- tree_ids[!(tree_ids %in% metadata_ids)]
  if (length(missing_in_metadata) > 0) {
    message <- paste(length(missing_in_metadata), "seqs missing from metadata:")
    stop(paste(message, paste0(head(missing_in_metadata, 5), collapse = ", ")))
  }
}

get_CI_date_min <- function(CI_date) {
  # tidytree::get_data returns node date ranges as list. Get the lower bound.
  if (is.null(CI_date)) {
    return(NA)
  }
  return(CI_date[1])
}

get_CI_date_max <- function(CI_date) {
  # tidytree::get_data returns node date ranges as list. Get the upper bound.
  if (is.null(CI_date)) {
    return(NA)
  }
  return(CI_date[2])
}

get_tips_under_node <- function(node, tree_data, n_tips) {
  # Base case: node is tip, return tip
  if (is_tip(node = node, n_tips = n_tips)) {
    return(node)
  }
  # Recursive case: paste together tips under children
  child_node_data <- get_child_node_data(node = node, tree_data = tree_data)
  return(unlist(lapply(FUN = get_tips_under_node, X = child_node_data$node, tree_data = tree_data, n_tips = n_tips)))
}

delete_internal_zero_branches <- function(tree_data, root_node, verbose = T) {
  # Delete zero-length internal branches 
  nodes_to_delete <- c()
  for (i in 1:nrow(tree_data)) {
    bl <- tree_data$branch.length[i]
    node_idx <- i
    while (bl == 0 & tree_data$node[node_idx] > root_node) {  # if a zero length branch b/w 2 internal nodes
      parent_idx <- which(tree_data$node == tree_data$parent[node_idx])
      bl <- tree_data[parent_idx, "branch.length"]
      node_idx <- parent_idx
      node <- tree_data[i, "node"]
      parent_node <- tree_data[parent_idx, "node"]
      if (verbose) {
        print(paste("Internal node", node, "is a daughter of a zero length branch"))
        print(paste("it will be deleted and its children will have parent node", parent_node))
      }
      tree_data[tree_data$parent == node, "parent"] <- parent_node
      nodes_to_delete <- c(nodes_to_delete, node)
    }
  }
  tree_data <- tree_data[!(tree_data$node %in% nodes_to_delete), ]
  return(tree_data)
}

get_locs_of_descendent_tips <- function(node, tree_data) {
  # Return unique locations amongst descendent tips of the node
  child_data <- get_child_node_data(node = node, tree_data = tree_data)
  return(unique(child_data$country_recoded[!is.na(child_data$country_recoded)]))
}

# Mapping from unambiguous canton codes to Nextstrain division names
# Compare https://github.com/nextstrain/ncov-ingest/blob/master/source-data/location_hierarchy.tsv
kanton_code_dict <- list(
  "AG" = "Aargau",
  "BE" = "Bern",
  "BL" = "Basel-Land",
  "BS" = "Basel-Stadt",
  "FR" = "Fribourg",
  "GE" = "Geneva",
  "GL" = "Glarus",
  "GR" = "Graubünden",
  "JU" = "Jura",
  "LU" = "Lucerne",
  "NW" = "Nidwalden",
  "OB" = "Obwalden",
  "SG" = "Sankt Gallen",
  "SH" = "Schaffhausen",
  "SO" = "Solothurn",
  "SZ" = "Schwyz",
  "TG" = "Thurgau",
  "TI" = "Ticino",
  "UR" = "Uri",
  "VD" = "Vaud",
  "VS" = "Valais",
  "ZG" = "Zug",
  "ZH" = "Zürich",
  # The following cantons are not in the Nextstrain dataset yet.
  # The mapping should be adapted once Nextstrain adds them.
  "AI" = "Appenzell Innerrhoden",
  "AR" = "Appenzell Ausserrhoden",
  "NE" = "Neuchâtel")

get_ns_division_from_canton_code <- function(canton_code) {
  if (!(canton_code %in% names(kanton_code_dict))) {
    stop(paste("Canton", canton_code, "not known."))
  }
  return(kanton_code_dict[canton_code])
}

get_CI_min <- function(CI) {
  CI_str <- gsub(
    x = gsub(x = CI, pattern = "c\\(", replacement = ""),
    pattern = "\\)", replacement = "")
  return(strsplit(CI_str, split = ", ")[[1]][1])
}
get_CI_max <- function(CI) {
  CI_str <- gsub(
    x = gsub(x = CI, pattern = "c\\(", replacement = ""),
    pattern = "\\)", replacement = "")
  return(strsplit(CI_str, split = ", ")[[1]][2])
}

get_week_since_epidemic_start <- function(date, first_seq_date = "2019-12-24") {
  # Return integer week for number of weeks passed since first_seq_date.
  first_seq_date <- as.Date(first_seq_date)
  days_since <- as.numeric(as.Date(date) - first_seq_date)
  weeks_since <- floor(days_since / 7)
  return(weeks_since)
}

get_days_since_epidemic_start <- function(date, first_seq_date = "2019-12-24") {
  # Return integer for number of days passed since first_seq_date.
  first_seq_date <- as.Date(first_seq_date)
  days_since <- as.numeric(as.Date(date) - first_seq_date)
  return(days_since)
}

get_date_range_for_week_since <- function(weeks_since, first_seq_date = "2019-12-24", date_format = "%b. %d", start_only = T) {
  # Return character string for start of week (or week date range) for integer weeks passed since first_seq_date.
  first_seq_date <- as.Date(first_seq_date)
  week_start <- first_seq_date + 7 * weeks_since 
  week_end <- week_start + 6
  if (start_only) {
    return(format(week_start, format = date_format))
  } else {
    return(paste(format(week_start, format = date_format), format(week_end, format = date_format), sep = " - "))
  }
}

get_weeks_since_to_date_range_mapping <- function(weeks_since, first_seq_date = "2019-12-24", date_format = "%b. %d", start_only = T) {
  # Return ordered data frame with one row per unique week in vector weeks_since, column giving start of week (or week date range) for that integer weeks passed since first_seq_date.
  first_seq_date <- as.Date(first_seq_date)
  complete_weeks_since <- min(weeks_since):max(weeks_since)
  complete_date_ranges <- get_date_range_for_week_since(
    weeks_since = complete_weeks_since, first_seq_date = first_seq_date, 
    date_format = date_format, start_only = start_only)
  return(data.frame(weeks_since = complete_weeks_since, date_range = complete_date_ranges))
}

get_chain_data_by_tip <- function(chains, metadata = NULL) {
  # Takes transmission chain data and spreads it into one row per tip with chain idx variable
  # Also merges in tip metadata
  chain_data_by_tip <- chains %>%
    mutate(chain_idx = 1:n()) %>%
    tidyr::separate(
      col = tips,
      sep = ", ",
      into = paste("tip", 1:max(chains$size), sep = "_"),
      remove = F,
      fill = "right") %>% 
    tidyr::pivot_longer(
      cols = paste("tip", 1:max(chains$size), sep = "_"),
      values_to = "tip",
      names_to = "tip_idx",
      names_prefix = "tip_") %>%
    filter(!(is.na(tip)))
  if (!is.null(metadata)) {
    chain_data_by_tip <- merge(
      x = chain_data_by_tip, y = metadata,
      by.x = "tip", by.y = "tree_label", all.x = T)
  }
  return(chain_data_by_tip)
}

#' Get number of sequences compared to confirmed cases per week.
#' @param qc_gisaid_query A database query providing sequence filtering criteria 
#' for seq data in table gisaid_sequence.
#' @param by_canton If True, returns results stratified by canton with additional column 'canton_code'.
#' @return Dataframe with one row per week and columns giving the number of 
#' sequences in the gisaid_sequence table from Switzerland passing the filter
#' criteria.
get_weekly_case_and_seq_data <- function(db_connection, qcd_gisaid_query, by_canton = F) {
  print("Getting weekly case and sequence data")
  sequence_data_query <- qcd_gisaid_query %>%
    filter(iso_country == "CHE") %>%
    mutate(
      is_viollier = originating_lab == "Viollier AG" & submitting_lab == "Department of Biosystems Science and Engineering, ETH Zürich",
      week = as.Date(date_trunc('week', date))) %>%
    group_by(is_viollier, week, division) %>%
    summarize(n_seqs = n(), .groups = "drop") %>%
    left_join(
      y = dplyr::tbl(db_connection, "swiss_canton") %>% 
        select(gisaid_division, canton_code),
      by = c("division" = "gisaid_division"))
  
  case_data_query <- dplyr::tbl(db_connection, "bag_test_numbers") %>%
    mutate(week = as.Date(date_trunc('week', date))) %>%
    group_by(week, canton) %>%
    summarise(n_conf_cases = sum(positive_tests, na.rm = T), .groups = "drop")
  
  weekly_case_and_seq_data <- dplyr::full_join(
    x = sequence_data_query, y = case_data_query, 
    by = c("week" = "week", "canton_code" = "canton")) %>% 
    collect() %>%
    tidyr::replace_na(replace = list(is_viollier = F)) %>%  # weeks with 0 seqs have all 0 seqs coming from other labs, dummy val
    tidyr::pivot_wider(
      names_from = "is_viollier", 
      values_from = "n_seqs", 
      names_prefix = "is_viollier_",
      values_fill = list(n_seqs = 0))
    if (!("is_viollier_TRUE" %in% colnames(weekly_case_and_seq_data))) {
      warning("There are no viollier sequences found in the Swiss dataset.")
      weekly_case_and_seq_data <- weekly_case_and_seq_data %>%
        mutate(is_viollier_TRUE = 0)
    }
    weekly_case_and_seq_data <- weekly_case_and_seq_data %>%
      rename("n_seqs_viollier" = "is_viollier_TRUE",
             "n_seqs_other" = "is_viollier_FALSE") %>%
      mutate(n_seqs_total = n_seqs_viollier + n_seqs_other,
             n_seqs_viollier = as.numeric(n_seqs_viollier),
             n_seqs_other = as.numeric(n_seqs_other),
             n_seqs_total = as.numeric(n_seqs_total),
             n_conf_cases = as.numeric(n_conf_cases)) %>%
      tidyr::replace_na(replace = list(n_conf_cases = 0,
                                       n_seqs_total = 0,
                                       n_seqs_other = 0))
  if (!by_canton) {
    weekly_case_and_seq_data <- weekly_case_and_seq_data %>%
      select(-c(canton_code, division)) %>%
      group_by(week) %>%
      summarise_all(sum) %>%
      mutate(canton_code = NA, division = NA)
  }
  return(weekly_case_and_seq_data)
}

#' Largest remainder method. 
#' Will take random from among possiblilities with equal reminders.
#' Also known as:  Hamilton, Hare-Niemeyer, Vinton method
#' Adapted from: https://github.com/polettif/proporz/blob/master/R/quota.R
#' @param votes Number of votes for each possibility
#' @param n_seats Number of seats to allocate to the possibilities
#' @return Number of seats allocated to each possibility
quota_largest_remainder = function(
  votes, n_seats, verbose = F
) { 
  if (sum(votes) == 0) {  # avoid undefined in quota when there are zero votes for any option
    seats_base = rep(0, length(votes))
    remainder = rep(0, length(votes))
    if (n_seats == 0) {
      return(rep(0, length(votes)))
    }
  } else {
    quota = n_seats * votes / sum(votes)
    seats_base = floor(quota)
    remainder = quota - seats_base
    if (all(remainder == 0)) {  # avoid 1:0 in calculating seats_rem evaluating to 1 and adding a spurious seat
      return(seats_base)
    }
  }
  n_seats_remaining = n_seats - sum(seats_base)
  seats_rem <- rep(0, length(votes))
  order_index = order(remainder, decreasing = TRUE)
  seats_rem[order_index[1:n_seats_remaining]] <- 1
  if (verbose) {
    print(paste("votes:", paste0(votes, collapse = ",")))
    print(paste("seats:", n_seats))
    print(paste("allocations:", paste0(seats_base + seats_rem, collapse = ",")))
  }
  return(seats_base + seats_rem)
}

#' Get symmetric limits for a plot
symmetric_limits <- function (x) {
  max <- max(abs(x))
  c(-max, max)
}

