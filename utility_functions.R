# Utility functions for project, mostly having to do with manipulating data and 
# trees.

#' Overwrites the open_database_connection function: Connection data (including the password) must be provided via
#' environment variables.
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
    password = connection_data$password,
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

get_n_descendent_focal_mrcas <- function(node, tree_data, chains, n_tips) {
  # Returns 1 if the node is a focal MRCA, otherwise returns the # focal MRCAs
  # descending from the node. Counts only the uppermost MRCA if there additional
  # focal MRCAs are nested under a focal MRCA higher in the tree.
  if (node %in% chains$ch_mrca) {
    return(1)
  } else if (node <= n_tips) {
    return(0)
  } else {
    child_node_data <- get_child_node_data(node = node, tree_data = tree_data)
  }
  n_descendent_swiss_mrcas <- sum(unlist(lapply(
    X = child_node_data$node,
    FUN = get_n_descendent_focal_mrcas,
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
#' @param subsample_by_canton If True, returns results stratified by canton with additional column 'canton_code'. Only available if focal country is CHE.
#' @param smooth_conf_cases If true, smooth confirmed case counts.
#' @return Dataframe with one row per week and columns giving the number of 
#' confirmed cases and sequences in the gisaid_sequence table from Switzerland 
#' passing the filter criteria, stratified into Viollier and non-Viollier samples.
#' Case-counts optionally smoothed across a 3-week window, reported at the cantonal level.
get_weekly_case_and_seq_data <- function(
  db_connection, qcd_gisaid_query, subsample_by_canton = F, smooth_conf_cases = F, focal_country
) {
  print("Getting weekly case and sequence data")
  weekly_case_and_seq_data <- query_weekly_case_and_seq_data(db_connection, qcd_gisaid_query, focal_country)
  
  print("Reformatting weekly case and sequence data")
  if (!("is_viollier_TRUE" %in% colnames(weekly_case_and_seq_data))) {
    warning("There are no viollier sequences found in the focal dataset. This is expected if focal country is not CHE.")
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
  if (!subsample_by_canton) {
    weekly_case_and_seq_data <- weekly_case_and_seq_data %>%
      select(-c(canton)) %>%
      group_by(week) %>%
      summarise_all(sum) %>%
      mutate(canton = focal_country)
  }
  
  if (smooth_conf_cases) {
    require(data.table)
    setDT(weekly_case_and_seq_data)     # converts test to a data.table in place
    setkey(weekly_case_and_seq_data,week, canton)
    weekly_case_and_seq_data <- with(weekly_case_and_seq_data, weekly_case_and_seq_data[order(canton, week),])
    weekly_case_and_seq_data[, n_conf_cases_smoothed := as.numeric(get.mav(n_conf_cases, 3)), by = canton]
    weekly_case_and_seq_data <- weekly_case_and_seq_data %>%
      mutate(n_conf_cases_raw = n_conf_cases,
             n_conf_cases = n_conf_cases_smoothed)
  }
  
  return(weekly_case_and_seq_data)
}

#' Thanks to: https://stackoverflow.com/questions/26198551/rolling-mean-moving-average-by-group-id-with-dplyr
#' @param bp A numeric vector.
#' @param n Window size.
#' @return Rolling average (left-aligned).
get.mav <- function(bp, n){
  require(zoo)
  bp <- na.locf(bp, na.rm = FALSE)  # replace NA with most recent non-NA prior to it
  if( length(bp) < n ) {
    return(bp)
  }
  # Keep first n - 1 values as-is
  return(c(bp[1:(n-1)], rollapply(bp, width = n, mean)))  
}

#' @return Dataframe with one row per week and columns giving the number of 
#' confirmed cases and sequences in the gisaid_sequence table from Switzerland 
#' passing the filter criteria, stratified into Viollier and non-Viollier samples.
#' If focal country is not CHE, no samples will be Viollier samples.
query_weekly_case_and_seq_data <- function(
  db_connection, qcd_gisaid_query, focal_country
) {
  sequence_data_query <- qcd_gisaid_query %>%
    filter(country == focal_country) %>%
    mutate(
      is_viollier = originating_lab == "Viollier AG" & submitting_lab == "Department of Biosystems Science and Engineering, ETH Zürich",
      week = as.Date(date_trunc('week', date))) %>%
    group_by(is_viollier, week, division) %>%
    summarize(n_seqs = n(), .groups = "drop")

  if (focal_country == "CHE") {
    sequence_data_query <- sequence_data_query  %>%
      left_join(
        y = dplyr::tbl(db_connection, "swiss_canton") %>%
          mutate(canton = canton_code) %>%
          select(gisaid_division, canton),
        by = c("division" = "gisaid_division"))
    case_data_query <- dplyr::tbl(db_connection, "bag_test_numbers") %>%
      mutate(week = as.Date(date_trunc('week', date))) %>%
      group_by(week, canton) %>%
      summarise(n_conf_cases = sum(positive_tests, na.rm = T), .groups = "drop")
  } else {
    case_data_query <- dplyr::tbl(db_connection, "ext_owid_global_cases") %>%
      filter(iso_country == focal_country) %>%
      mutate(week = as.Date(date_trunc('week', date))) %>%
      group_by(week) %>%
      summarise(n_conf_cases = sum(new_cases, na.rm = T), .groups = "drop") %>%
      mutate(canton = NA)
  }
  
  weekly_case_and_seq_data <- dplyr::full_join(
    x = sequence_data_query %>%
      tidyr::replace_na(replace = list(canton = focal_country)),
    y = case_data_query %>%
      tidyr::replace_na(replace = list(canton = focal_country))) %>%  # can't join NA with NA, need to replace NA with country code
    collect() %>%
    tidyr::replace_na(replace = list(is_viollier = F)) %>%  # weeks with 0 seqs have all 0 seqs coming from other labs, dummy val
    tidyr::pivot_wider(
      names_from = "is_viollier", 
      values_from = "n_seqs", 
      names_prefix = "is_viollier_",
      values_fill = list(n_seqs = 0))

  # Lump sequences with division not mapping to a canton (e.g. 'Basle') into NA canton for that week
  weekly_case_and_seq_data <- weekly_case_and_seq_data %>%
    select(-division) %>%
    group_by(week, canton) %>%
    summarise_all(sum)
  
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


# Thanks to https://github.com/YuLab-SMU/treeio/blob/master/R/beast.R
# Need to use this version of the function because in-package version doesn't handle the ambigous "2021" date output by LSD well
read.stats_beast_internal <- function(beast, tree) {
  ##tree <- gsub(" ", "", tree)
  ## tree2 <- gsub("\\[[^\\[]*\\]", "", tree)
  ## phylo <- read.tree(text = tree2)
  ## tree2 <- add_pseudo_nodelabel(phylo, tree2)
  
  phylo <- read.tree(text = tree)
  tree2 <- add_pseudo_nodelabel(phylo)
  
  ## node name corresponding to stats
  nn <- strsplit(tree2, split=",") %>% unlist %>%
    strsplit(., split="\\)") %>% unlist %>%
    gsub("\\(*", "", .) %>%
    gsub("[:;].*", "", .) %>%
    gsub(" ", "", .) %>%
    gsub("'", "", .) %>%
    gsub('"', "", .)
  
  phylo <- read.tree(text = tree2)
  root <- rootnode(phylo)
  nnode <- phylo$Nnode
  
  tree_label <- c(phylo$tip.label, phylo$node.label)
  ii <- match(nn, tree_label)
  
  if (any(grepl("TRANSLATE", beast, ignore.case = TRUE))) {
    label2 <- c(phylo$tip.label,
                root:getNodeNum(phylo))
    ## label2 <- c(treeinfo[treeinfo$isTip, "label"],
    ##             root:(root+nnode-1))
    
  } else {
    ## node <- as.character(treeinfo$node[match(nn, treeinfo$label)])
    label2 <- as.character(1:getNodeNum(phylo))
  }
  node <- label2[match(nn, tree_label)]
  
  ## stats <- unlist(strsplit(tree, "\\["))[-1]
  ## stats <- sub(":.+$", "", stats
  ## BEAST1 edge stat fix
  tree <- gsub("\\]:\\[&(.+?\\])", ",\\1:", tree)
  # t1:[&mutation="test1"]0.04 -> t1[&mutation="test1"]:0.04
  tree <- gsub(":(\\[.+?\\])", "\\1:", tree)
  
  if (grepl("\\:[0-9\\.eEL+\\-]*\\[", tree) || grepl("\\]\\[", tree)){
    # t1:0.04[&mutation="test1"] -> t1[&mutation="test1"]:0.04
    # or t1[&prob=100]:0.04[&mutation="test"] -> t1[&prob=100][&mutation="test"]:0.04 (MrBayes output)
    # pattern <- "(\\w+)?(:?\\d*\\.?\\d*[Ee]?[\\+\\-]?\\d*)?(\\[&.*?\\])"
    pattern <- "(\\w+)?(:\\d*\\.?\\d*[Ee]?[\\+\\-]?\\L*\\d*)?(\\[&.*?\\])"
    tree <- gsub(pattern, "\\1\\3\\2", tree)
  }
  #if (grepl("\\]:[0-9\\.eE+\\-]*\\[", tree) || grepl("\\]\\[", tree)) {
  #    ## MrBayes output
  #    stats <- strsplit(tree, "\\]:[0-9\\.eE+\\-]*\\[") %>% unlist
  #    lstats <- lapply(stats, function(x) {
  #        unlist(strsplit(x, split="\\][,\\)]"))
  #    })
  #    
  #    for (i in seq_along(stats)) {
  #        n <- length(lstats[[i]])
  #        if (i == length(stats)) {
  #            stats[i] <- lstats[[i]][n]
  #        } else {
  #            stats[i] <- paste0(lstats[[i]][n],
  #                               sub("&", ",", lstats[[i+1]][1])
  #            )
  #        }
  #    }
  #    stats <- gsub("\\]\\[&", ",", stats)
  #} else {
  #    ## BEAST output
  #    stats <- strsplit(tree, ":") %>% unlist
  #}
  stats <- strsplit(tree, ":") %>% unlist
  names(stats) <- node
  
  stats <- stats[grep("\\[", stats)]
  stats <- sub("[^\\[]*\\[", "", stats)
  
  stats <- sub("^&", "", stats)
  # this is for MrBayes output 
  stats <- sub("\\]\\[&", ",", stats)
  stats <- sub("];*$", "", stats)
  stats <- gsub("\"", "", stats)
  
  stats2 <- lapply(seq_along(stats), function(i) {
    x <- stats[[i]]
    y <- unlist(strsplit(x, ","))
    # the stats information does not has always {}
    #sidx <- grep("=\\{", y)
    #eidx <- grep("\\}$", y)
    # [&mutation="test1,test2",rate=80,90]
    sidx1 <- grep("=", y)
    eidx1 <- sidx1 - 1
    eidx1 <- c(eidx1[-1], length(y))
    # for better parsing [&mutation="test",name="A"] single value to key.
    sidx <- sidx1[!(sidx1==eidx1)]
    eidx <- eidx1[!(sidx1==eidx1)]
    
    flag <- FALSE
    if (length(sidx) > 0) {
      flag <- TRUE
      SETS <- lapply(seq_along(sidx), function(k) {
        p <- y[sidx[k]:eidx[k]]
        gsub(".*=\\{", "", p) %>% 
          gsub("\\}$", "", .) %>%
          gsub(".*=", "", .)
      })
      names(SETS) <- gsub("=.*", "", y[sidx])
      
      kk <- lapply(seq_along(sidx), function(k) {
        sidx[k]:eidx[k]
      }) %>%
        unlist
      y <- y[-kk]
    }
    
    if (length(y) == 0)
      return(SETS)
    
    name <- gsub("=.*", "", y)
    val <- gsub(".*=", "", y) %>%
      gsub("^\\{", "", .) %>%
      gsub("\\}$", "", .)
    
    if (flag) {
      nn <- c(name, names(SETS))
    } else {
      nn <- name
    }
    
    res <- rep(NA, length(nn))
    names(res) <- nn
    
    for (i in seq_along(name)) {
      res[i] <- if(is_numeric(val[i])) as.numeric(val[i]) else val[i]
    }
    if (flag) {
      j <- i
      for (i in seq_along(SETS)) {
        if(is_numeric(SETS[[i]])) {
          res[i+j] <- list(as.numeric(SETS[[i]]))
        } else {
          res[i+j] <- SETS[i]
        }
      }
    }
    
    return(res)
  })
  
  nn <- lapply(stats2, names) %>% unlist %>%
    unique %>% sort
  
  
  stats2 <- lapply(stats2, function(x) {
    y <- x[nn]
    names(y) <- nn
    y[vapply(y, is.null, logical(1))] <- NA
    y
  })
  
  stats3 <- do.call(rbind, stats2)
  stats3 <- as_tibble(stats3)
  
  ## no need to extract sd from prob+-sd
  ## as the sd is stored in prob_stddev
  ##
  ## "prob_stddev"   "prob(percent)" "prob+-sd"
  ##
  ##
  ##
  ## idx <- grep("\\+-", colnames(stats3))
  ## if (length(idx)) {
  ##     for (i in idx) {
  ##         stats3[,i] <- as.numeric(gsub("\\d+\\+-", "", stats3[,i]))
  ##     }
  ## }
  
  cn <- gsub("(\\d+)%", "0.\\1", colnames(stats3))
  cn <- gsub("\\(([^\\)]+)\\)", "_\\1", cn)
  ## cn <- gsub("\\+-", "_", cn)
  
  colnames(stats3) <- cn
  stats3$node <- names(stats)
  
  i <- vapply(stats3,
              function(x) max(vapply(x, length, numeric(1))),
              numeric(1))
  
  for (j in which(i==1)) {
    stats3[,j] <- unlist(stats3[,j])
  }
  stats3$node <- as.integer(stats3$node)
  return(stats3)
}


add_pseudo_nodelabel <- function(phylo) {
  if(is.null(phylo$node.label)) {
    nnode <- phylo$Nnode
    phylo$node.label <- paste("X", 1:nnode, sep="")
    ## for (i in 1:nnode) {
    ##     treetext <- sub("\\)([:;])",
    ##                     paste0("\\)", nlab[i], "\\1"),
    ##                     treetext)
    ## }
  }
  ## if tip.label contains () which will broken node name extraction
  phylo$tip.label <- gsub("[\\(\\)]", "_", phylo$tip.label)
  
  treetext <- write.tree(phylo)
  return(treetext)
}
