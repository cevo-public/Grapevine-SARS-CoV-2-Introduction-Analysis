require(dplyr)

get_week_to_date_range_mapping <- function(min_date, max_date, date_format = "%b. %d") {
  week_to_date_ranges <- data.frame(
    # Generate mapping from weeks to date ranges
    date = seq(from = as.Date(min_date), to = as.Date(max_date), by = 1)) %>%
    mutate(week = week(date)) %>%
    group_by(week) %>%
    summarize(date_range = paste(
      format(min(date), date_format), 
      format(max(date), date_format), sep = " - "))
  return(week_to_date_ranges)
}

kanton_name_dict <- list(
  "ZH" = "Zurich",
  "ZH" = "Zürich",
  "BE" = "Bern",
  "LU" = "Luzern",
  "LU" = "Lucerne",
  "UR" = "Uri",
  "SZ" = "Schwyz",
  "OW" = "Obwalden",
  "NW" = "Nidwalden",
  "GL" = "Glarus",
  "ZG" = "Zug",
  "FR" = "Fribourg",
  "SO" = "Solothurn",
  "BS" = "Basel-Stadt",
  "BL" = "Basel-Landschaft",
  "BL" = "Basel-Land",
  "SH" = "Schaffhausen",
  "AR" = "Appenzell Ausserrhoden",
  "AI" = "Appenzell Innerrhoden",
  "SG" = "St Gallen",
  "SG" = "Sankt Gallen",
  "GR" = "Graubünden",
  "AG" = "Aargau",
  "TG" = "Thurgau",
  "TI" = "Ticino",
  "VD" = "Vaud",
  "VS" = "Valais",
  "NE" = "Neuchâtel",
  "GE" = "Geneva",
  "JU" = "Jura")

get_kanton_code_from_fullnames <- function(canton_fullname) {
  if (!(canton_fullname %in% kanton_name_dict)) {
    stop(paste("Canton", canton_fullname, "not known."))
  }
  return(names(kanton_name_dict)[which(kanton_name_dict == canton_fullname)])
}

load_metadata_files <- function(path, pattern) {
  metadata_files <- list.files(
    path = path, 
    pattern = pattern, 
    full.names = T, 
    recursive = T)
  is_first <- T
  for (file in metadata_files) {
    metadata_temp <- read.delim(file = file, stringsAsFactors = F)
    prefix = strsplit(file, split = "/")[[1]][10]
    metadata_temp$`Tree size` <- strsplit(prefix, split = "_")[[1]][2]
    metadata_temp$Replicate <- strsplit(prefix, split = "_")[[1]][4]
    if (is_first) {
      is_first <- F
      metadata_with_run <- metadata_temp
    } else {
      metadata_with_run <- rbind(metadata_with_run, metadata_temp)
    }
  }
  return(metadata_with_run)
}

load_cluster_data <- function(path, pattern) {
  files <- list.files(path = path, pattern = pattern, full.names = F)
  is_first <- T
  for (file in files) {
    cluster_data_temp <- read.delim(
      file = paste(path, file, sep = "/"), stringsAsFactors = F)
    prefix = strsplit(file, split = "_clusters.txt")[[1]][1]
    cluster_data_temp$run <- prefix
    if (is_first) {
      is_first <- F
      cluster_data <- cluster_data_temp
    } else {
      cluster_data <- rbind(cluster_data, cluster_data_temp)
    }
  }
  date_col_filter <- grepl(x =  colnames(cluster_data), pattern = "_tmrca")
  date_cols <- colnames(cluster_data)[date_col_filter]
  cat("interpreting the following as date columns:\n")
  cat(paste0(date_cols, collapse = "\n"))
  for (col in date_cols) {
    cluster_data[[col]] <- as.Date(cluster_data[[col]])
  }
  return(cluster_data)
}

extract_small_clade <- function(tree, min_size = 5, max_size = 10, max_tries = 1000, ...) {
  # Choose random clades on tree until one of the appropriate size is found
  tree_size <- max_size + 1
  n_tries <- 1
  while((tree_size > max_size | tree_size < min_size) & n_tries < max_tries) {
    node <- sample(1:length())
    tree_small <- ape::extract.clade(phy = tree, node = sample(tree$node.label, size = 1))
    tree_size <- length(tree_small$tip.label)
    n_tries <- n_tries + 1
  }
  if (n_tries < max_tries) {
    return(tree_small)
  } else {
    stop("Timed out trying to find a clade of the correct size.")
  }
}

extract_outgroup_and_random_tips <- function(tree, size = 10, outgroup_tips) {
  # Choose outgroup seqs and a selection of random other tips from tree
  random_tips <- sample(x = tree$tip.label, size = size - 2, replace = F)
  tree_small <- ape::keep.tip(phy = tree, tip = c(outgroup_tips, random_tips))
  return(tree_small)
}

extract_clade_by_node_name <- function(tree, node = 3360, ...) {
  # Extract node from tree
  # tree_small <- treeio::tree_subset(tree = tree, node = node, levels_back = 0)
  tree_small <- ape::extract.clade(phy = tree, node = node)
  return(tree_small)
}

extract_outgroup_clade <- function(tree, outgroup_tips, ...) {
  # Extract clade containing outgroup seqs
  mrca <- ape::getMRCA(phy = tree, tip = outgroup_tips)
  tree_small <- ape::extract.clade(phy = tree, node = mrca)
  return(tree_small)
}

get_date_from_tip_label <- function(tip_label) {
  return(strsplit(tip_label, split = "\\|")[[1]][3])
}
get_loc_from_tip_label <- function(tip_label) {
  return(strsplit(
    strsplit(tip_label, split = "\\|")[[1]][1],
    split = "_")[[1]][1])
}

get_cluster_data_by_tip <- function(cluster_data, metadata) {
  # Takes cluster data and spreads it into one row per tip with cluster idx variable
  # Also merges in tip metadata
  cluster_data_by_tip <- cluster_data %>%
    mutate(cluster_idx = 1:n()) %>%
    tidyr::separate(
      col = tips,
      sep = ", ",
      into = paste("tip", 1:max(cluster_data$size), sep = "_"),
      remove = F,
      fill = "right") %>% 
    tidyr::pivot_longer(
      cols = paste("tip", 1:max(cluster_data$size), sep = "_"),
      values_to = "tip",
      names_to = "tip_idx",
      names_prefix = "tip_") %>%
    filter(!(is.na(tip)))
  cluster_data_by_tip <- merge(
    x = cluster_data_by_tip, y = metadata,
    by.x = "tip", by.y = "strain", all.x = T)
  return(cluster_data_by_tip)
}

# Create layout where cantons are placed according to geographic location
canton_loc_data <- data.frame(
  cantons = c(
    "Aargau", "Basel-Land", "Basel-Stadt", "Bern", "Fribourg", "Geneva", 
    "Glarus", "Graubünden", "Jura", "Lucerne", "Sankt Gallen", "Schaffhausen",
    "Schwyz", "Solothurn", "Thurgau", "Ticino", "Uri", "Valais", "Vaud", "Zug", 
    "Zürich"),
  y = c(47.3909747, 47.4831796, 47.5545913, 46.9546485, 46.8031637, 46.2050241, 47.0132885, 46.61563, 47.3271569, 47.0547335, 47.4240435, 47.7152496,
        47.0257335,47.2080787,47.5352122,46.2237891,46.7600964,46.2547441,46.585495,47.1517046,47.3774336),
  x = c(8.0347348, 7.6919098, 7.5593547, 7.3246586, 7.1422125, 6.1089833, 8.9150336, 9.0098434, 6.9184466, 8.2120936, 9.2931586, 8.6107398,
        8.6208786, 7.5131775, 8.7913418, 8.2091573, 8.3968013, 7.0627024, 6.0948383, 8.4810097, 8.4665036))

shift_coordinate <- function(x, start_old, end_old, start_new, end_new) {
  return(start_new + ((end_new - start_new) / (end_old - start_old)) * (x - start_old))
}

canton_loc_data$x_new <- unlist(lapply(
  FUN = shift_coordinate, 
  X = canton_loc_data$x,
  start_old = min(canton_loc_data$x),
  end_old = max(canton_loc_data$x),
  start_new = -1, end_new = 1))
canton_loc_data$y_new <- unlist(lapply(
  FUN = shift_coordinate, 
  X = canton_loc_data$y,
  start_old = min(canton_loc_data$y),
  end_old = max(canton_loc_data$y),
  start_new = -1, end_new = 1))

# Shift Basel-Stadt away from Basel-Land a bit so the arrow is visible
canton_loc_data[3, "x_new"] <- canton_loc_data[3, "x_new"] - 0.08
canton_loc_data[3, "y_new"] <- canton_loc_data[3, "y_new"] + 0.08

# Sequences that were published by us (author is "Christian Biesel et al") but 
# have non-ETH & Violleri authors b/c they come from UHZ
UHZ_SAMPLES <- c(
  "EPI_ISL_483668", "EPI_ISL_483669", "EPI_ISL_483670", "EPI_ISL_483671", 
  "EPI_ISL_483672", "EPI_ISL_483673", "EPI_ISL_483674", "EPI_ISL_483675", "EPI_ISL_483676", 
  "EPI_ISL_483677", "EPI_ISL_483678", "EPI_ISL_483679", "EPI_ISL_483680", "EPI_ISL_483681", 
  "EPI_ISL_483682", "EPI_ISL_483683", "EPI_ISL_483684", "EPI_ISL_483685")

loadLogFiles <- function(logFiles, datesFile, sampling, sequences, burninFrac=0.1) {
  require(tidytable)
  require(dplyr)
  require(tidyverse)
  # Load BDKSY logfiles, output list of R0 values and sampling proportion values
  # Written by Tim!
  
  dates <- read.delim(datesFile)$date
  mostRecentSample <- max(as.Date(dates))
  
  data <- NULL
  for (file in logFiles) {
    cat(paste("Loading", file, "..."))
    df <- read.delim(file)
    N <- dim(df)[1]
    data <- bind_rows(data, df[-(1:ceiling(N*burninFrac)),])
    cat(" done.\n")
  }
  
  renamer <- function(colname) {
    indices <- as.numeric(str_split(colname, "\\.", simplify=TRUE)[,-1])
    output <- col_dates[indices+1]
    return(output)
  }
  col_dates <- as.character(rev(seq.Date(mostRecentSample-lubridate::weeks(30),
                                         mostRecentSample, by="week")))
  cols <- paste0("R0Values.", 0:30)
  
  R0Data <- data %>%
    select("Sample", cols) %>%
    tidytable::rename_with.(renamer, starts_with("R0Values.")) %>%
    pivot_longer(all_of(col_dates), names_to="Date", values_to="R0") %>%
    mutate(Date=lubridate::ymd(Date)) %>%
    mutate(Sampling=sampling, Sequences=sequences)
  
  
  if (sampling == "Skyline") {
    cols <- paste0("sampValues.", 0:30)
    samplingPropData <- data %>%
      select("Sample", cols) %>%
      tidytable::rename_with.(renamer, starts_with("sampValues.")) %>%
      pivot_longer(all_of(col_dates), names_to="Date", values_to="sampProp") %>%
      mutate(Date=lubridate::ymd(Date)) %>%
      mutate(Sampling=sampling, Sequences=sequences)
  } else {
    samplingPropData <- NULL
  }
  
  return(list(R0Data=R0Data, samplingPropData=samplingPropData))
}
