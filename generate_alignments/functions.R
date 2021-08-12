#' Direct sample selection to appropriate function.
#' @param focal_country ISO code for focal country.
select_sequences <- function(
  qcd_gisaid_query, db_connection, max_sampling_frac, focal_country = 'CHE', 
  verbose = T, subsample_by_canton = T, smooth_conf_cases = F,
  outdir = NULL, max_date, min_date
) {
    all_samples <- qcd_gisaid_query %>%
        filter(country == focal_country) %>%
        mutate(week = as.Date(date_trunc('week', date))) %>%
        select(strain, week, division, gisaid_epi_isl) %>%
        left_join(
          y = dplyr::tbl(db_connection, "swiss_canton") %>%
            mutate(canton = canton_code) %>%
            select(gisaid_division, canton),
          by = c("division" = "gisaid_division")) %>%
        collect()
    print(paste("There are a total of", nrow(all_samples), "focal samples passing the qcd_gisaid_query"))
    if (max_sampling_frac > 0 & max_sampling_frac <= 1) {
      downsample_data <- downsample_focal_sequences(
        qcd_gisaid_query = qcd_gisaid_query,
        all_samples = all_samples,
        db_connection = db_connection,
        max_sampling_frac = max_sampling_frac,
        focal_country = focal_country,
        verbose = verbose,
        subsample_by_canton = subsample_by_canton,
        smooth_conf_cases = smooth_conf_cases,
        outdir = outdir,
        max_date = max_date,
        min_date = min_date)
      return(downsample_data)
    } else if (max_sampling_frac == -1) {
      print("Not downsampling sequences.")
      return(qcd_gisaid_query)
    } else if (max_sampling_frac < 0 | max_sampling_frac > 1) {
      stop("Error in selecting sequences: max_sampling_frac outside of range 0-1 or -1 for all samples.")
    } else {
      stop("Error in selecting sequences: unspecified input error.")
    }
}

#' Downsample Swiss sequences on GISAID proportionally to weekly confirmed cases.
#' Optionally downsample proportionally to weekly confirmed cases in each canton.
#' If there aren't enough sequences in a week from the region, all sequences are taken.
#' When confirmed cases are not attributed to a region, the proportional number of sequences are taken without regard to region.
#' Optionally favor samples with recorded foreign exposure information per BAG meldeformular or GISAID metadata before filling up the week's sequence quota with non-exposed/unknown samples.
#' @param qcd_gisaid_query Query of database table gisaid_sequence with QC filters.
#' @param max_sampling_frac Take no more than this fraction of confirmed cases each week.
#' @param subsample_by_canton If true, take sequences proportional to confirmed cases per canton. Only available if focal country is CHE.
#' @param smooth_conf_cases If true, use smoothed confirmed case counts.
#' @param output_summary_table If true, output table summarizing sampled sequences per region per week compare to confirmed cases
#' @return qcd_gisaid_query Query filtered to include foreign samples & only selected swiss samples
downsample_focal_sequences <- function (
  qcd_gisaid_query, all_samples, db_connection, max_sampling_frac,
  verbose = T, subsample_by_canton = T, smooth_conf_cases = F, outdir = NULL, focal_country, min_date, max_date
) {
  print("Downsampling focal sequences")

  # Get available seqs and cases numbers per week, possibly stratified by canton, possibly smoothed
  sampling_data_raw <- get_weekly_case_and_seq_data(
    db_connection = db_connection,
    qcd_gisaid_query = qcd_gisaid_query,
    subsample_by_canton = subsample_by_canton,
    focal_country = focal_country,
    smooth_conf_cases = smooth_conf_cases) %>%
    filter(week <= max_date, week >= min_date)

  # Calculate # seqs per week
  # and divide proportionally amongst cantons if data stratified by canton
  sampling_data <- sampling_data_raw %>%
    filter(canton != "FL") %>%  # don't count samples from Liechtenstein
    group_by(week) %>%
    mutate(
      max_seqs_from_week = floor(sum(n_conf_cases) * max_sampling_frac),
      n_ideal_sample = quota_largest_remainder(
        votes = n_conf_cases,
        n_seats = max_seqs_from_week[1]),
      n_to_sample = case_when(
        canton == focal_country ~ 0,  # placeholder, to be filled below
        T ~ pmin(n_ideal_sample, n_seqs_total)))

  # For weeks where some case data not aggregated by canton, take leftover sequences from whole-country
  leftover_seqs_by_week <- sampling_data %>%
    group_by(week) %>%
    summarize(n_leftover_seqs = sum(n_seqs_total) - sum(n_to_sample)) %>%  # includes extra sequences from cantons and seqs not assigned a canton
    mutate(canton = focal_country)
  sampling_data <- merge(x = sampling_data, y = leftover_seqs_by_week, all.x = T, by = c("week", "canton"))
  sampling_data <- sampling_data %>%
    mutate(
      n_to_sample = case_when(
        canton == focal_country ~ pmin(n_ideal_sample, n_leftover_seqs),
        T ~ n_to_sample))

  # select sequences from cantons first, so that non-canton-specific sampling doesn't take too many first
  canton_levels <- c(sort(unique(sampling_data$canton[sampling_data$canton != focal_country])), focal_country)
  sampling_data$canton<- factor(
    x = sampling_data$canton,
    levels = canton_levels)
  sampling_data <- sampling_data %>% arrange(canton, week)

  # sample sequences
  sampled_strains <- c()
  for (i in 1:nrow(sampling_data)) {
    week_i <- sampling_data[[i, "week"]]
    canton_i <- as.character(sampling_data[[i, "canton"]])
    n_samples_i <- as.numeric(sampling_data[[i, "n_to_sample"]])
    
    if (n_samples_i == 0) {  # if no samples desired, continue
      next
    }
    
    if (canton_i == focal_country) {  # if not stratifying by canton, or the canton is not known, take randomly from all Switzerland
      all_samples_i <- all_samples %>%
        filter(week == week_i, !(strain %in% !! sampled_strains))
    } else {  # else take from the right canton
      all_samples_i <- all_samples %>%
        filter(week == week_i, canton == canton_i, !(strain %in% !! sampled_strains))
    }
    if (verbose) {
      print(paste(
        "sampling", n_samples_i, 
        "samples out of", nrow(all_samples_i), 
        "available from", ifelse(test = is.na(canton_i), yes = paste("all", focal_country), no = canton_i),
        "in week", format(as.Date(week_i), "%Y-%m-%d")))
    }
    sampled_strains_i <- sample_strains(all_samples_i, n_samples_i)
    sampled_strains <- c(sampled_strains, sampled_strains_i)
  }

  # output downsampling data table and figure
  if (!(is.null(outdir))) {
    report_downsampling(sampling_data, outdir, max_sampling_frac, focal_country)
  }
  
  # update master seq query to exclude non-sampled focal strains
  qcd_gisaid_query <- qcd_gisaid_query %>%
    filter(country != !! focal_country | strain %in% !! sampled_strains)

  return(qcd_gisaid_query)
}

sample_strains <- function(all_samples_i, n_samples_i) {
  shuffled_samples_i <- all_samples_i[sample(nrow(all_samples_i), replace = F), ]
  sampled_strains_i <- shuffled_samples_i$strain[1:n_samples_i]
  return(sampled_strains_i)
}

report_downsampling <- function(sampling_data, outdir, max_sampling_frac, focal_country) {
  write.csv(
    x = sampling_data, 
    file = paste(outdir, paste0(focal_country, "_downsampling_data.csv"), sep = "/"),
    row.names = F)
  sampling_plot <- ggplot(
    data = sampling_data %>% mutate(canton = tidyr::replace_na(data = canton, replace = focal_country)),
    aes(x = as.Date(week))) + 
    geom_col(aes(y = n_to_sample, fill = "Number of sequences analyzed")) + 
    geom_col(aes(y = -n_conf_cases * max_sampling_frac, 
                 fill = paste(max_sampling_frac * 100, "% of confirmed cases", sep = ""))) + 
    facet_wrap(canton ~ ., scales = "free_y") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b. %y") +
    scale_y_continuous(limits = symmetric_limits) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_blank()) + 
    labs(x = element_blank(), y = "Count")
  ggsave(
    filename = paste(outdir, paste0(focal_country, "_downsampling.png"), sep = "/"),
    plot = sampling_plot)
}

# ------------------------------------------------------------------------------

#' QC gisaid data and split by pangolin lineage, aggregating lineages that are 
#' predominantly focal into the parent lineage.
get_pangolin_lineages <- function(db_connection, outdir, qcd_gisaid_query, focal_country) {
  
  print("Querying database for pangolin lineages.")
  pangolin_lineages <- qcd_gisaid_query %>%
    mutate(is_focal = country == focal_country) %>%
    group_by(pangolin_lineage, is_focal) %>%
    summarize(n_seqs = n(), .groups = "drop") %>%
    collect()

  print("Querying database for pangolin alias lookup table.")
  alias_table <- dplyr::tbl(db_connection, "pangolin_lineage_alias") %>% collect()

  print("Grouping by expanded pangolin lineage names, in case any are aliased.")
  pangolin_lineages <- pangolin_lineages %>%
    mutate(pangolin_lineage = unlist(lapply(
      FUN = expand_pangolin_lineage_aliases,
      X = pangolin_lineage,
      alias_table = alias_table))) %>%
    group_by(pangolin_lineage, is_focal) %>%
    summarize(n_seqs = sum(n_seqs), .groups = "drop")

  print("Storing alias information for previously aliased lineages.")
  pangolin_lineages <- pangolin_lineages %>% 
    mutate(
      lineage_alias = unlist(lapply(
        FUN = get_pangolin_lineage_aliases,
        X = pangolin_lineage,
        alias_table = alias_table)))
  
  print("Aggregating predominantly focal lineages into parent lineage.")
  pangolin_lineages_aggregated <- pangolin_lineages %>%
    tidyr::pivot_wider(
      names_from = is_focal, names_prefix = "is_focal_", values_from = n_seqs) %>%
    tidyr::replace_na(replace = list("is_focal_TRUE" = 0, "is_focal_FALSE" = 0)) %>%  # because by default 0s show up as NA after group_by %>% summarize
    filter(is_focal_TRUE > 0) %>%
    mutate(
      should_aggregate = is_focal_TRUE > is_focal_FALSE,
      parent_lineage = unlist(lapply(
        FUN = get_parent_lineage,
        X = pangolin_lineage)),
      n_lineages_aggregated = 1,
      lineages_aggregated = case_when(
        is.na(lineage_alias) ~ pangolin_lineage,
        T ~ paste(pangolin_lineage, lineage_alias, sep = ", ")))

  while(any(!is.na(pangolin_lineages_aggregated$parent_lineage) & 
            pangolin_lineages_aggregated$should_aggregate)) {
    pangolin_lineages_aggregated <- aggregate_predominatly_focal_lineages(
      pangolin_lineages_aggregated, alias_table)
  }
  
  write.csv(
    x = pangolin_lineages_aggregated,
    file = paste(outdir, "pangolin_lineages_aggregated.csv", sep = "/"),
    row.names = F)
  return(pangolin_lineages_aggregated)
}

#' @return Expansion of pangolin lineage alias, if alised
expand_pangolin_lineage_aliases <- function(pangolin_lineage, alias_table) {
  prefix <- stringr::str_extract(string = pangolin_lineage, pattern = "^[[:alpha:]]*")
  suffix <- gsub(x = pangolin_lineage, pattern = "^[[:alpha:]]*", replacement = "")
  if (prefix %in% alias_table$alias) {
    new_prefix <- unlist(unname(alias_table[alias_table$alias == prefix, "full_name"]))
  } else {
    new_prefix <- prefix  # not an alias, return same lineage
  }
  return(paste0(new_prefix, suffix))                   
}

#' @return Contraction of pangolin lineage into an alias, if any alias exists
#' E.g. given pangolin_lineage 'B.1.177.15.1', return 'AA.1' because 'B.1.177.15' = 'AA'
get_pangolin_lineage_aliases <- function(pangolin_lineage, alias_table, verbose = F) {
  suffix <- ""
  while (!is.na(pangolin_lineage)) {
    if (verbose) print(paste("Trying prefix", pangolin_lineage, "with suffix", suffix))
    if (pangolin_lineage %in% alias_table$full_name) {
      new_prefix <- unlist(unname(alias_table[alias_table$full_name == pangolin_lineage, "alias"]))
      if (suffix != "") {
        alias <- paste(new_prefix, suffix, sep = ".")
      } else {
        alias <- new_prefix
      }
      return(alias)
    }
    suffix_addition <- sub("^(.*)\\.", replacement = "", pangolin_lineage)
    if (suffix != "") {
      suffix <- paste(suffix_addition, suffix, sep = ".")
    } else {
      suffix <- suffix_addition
    }
    pangolin_lineage <- tryCatch(
    {get_parent_lineage(pangolin_lineage)},
    warning = function(cond) {return(NA)})  # catch warning that there's no more valid parent lineage silently
  }
  return(NA)
}

#' @param pangolin_lineage Lineage name, assumes alias are already expanded.
#' @return parent lineage, or NA if there is no valid parent lineage.
get_parent_lineage <- function(pangolin_lineage) {
  if (grepl(x = pangolin_lineage, pattern = "^X")) {  
    # Recombinants don't have parents
    print(paste("No valid parent lineage for recombinant", pangolin_lineage, "\n"))
    return(NA)
  } else if (grepl(x = pangolin_lineage, pattern = "\\.[[:digit:]]*\\.")) {  
    # Given B.1.1.1, return B.1.1; Given B.1.1 return B.1
    return(sub(".[^.]+$", "", pangolin_lineage))
  } else if (grepl(x = pangolin_lineage, pattern = "(A|B)\\.[[:digit:]]*$")) {
    # Special case, only A and B are valid single-letter lineages
    return(sub(".[^.]+$", "", pangolin_lineage))
  } else {
    warning(paste("Cannot find a valid parent lineage for", pangolin_lineage, "\n"))
    return(NA)
  }
}

#' For pangolin lineages where > 50% of samples are focal, assign samples to the
#' parent lineage. The goal here is that the MRCA of all lineages is outside of 
#' focal country.
aggregate_predominatly_focal_lineages <- function(pangolin_lineages, alias_table) {
  pangolin_lineages_aggregated <- pangolin_lineages %>% 
    mutate(
      pangolin_lineage_rollup = case_when(
        is_focal_TRUE > is_focal_FALSE & !is.na(parent_lineage) ~ parent_lineage,
        T ~ pangolin_lineage),
      is_high_level_focal_lineage = should_aggregate & is.na(parent_lineage)) %>%
    group_by(pangolin_lineage_rollup) %>%
    summarize(is_focal_TRUE = sum(is_focal_TRUE),
              is_focal_FALSE = sum(is_focal_FALSE),
              n_lineages_aggregated = sum(n_lineages_aggregated),
              lineages_aggregated = paste0(lineages_aggregated, collapse = ", "))
  
  pangolin_lineages_aggregated <- pangolin_lineages_aggregated %>%
    mutate(pangolin_lineage = pangolin_lineage_rollup,
           should_aggregate = is_focal_TRUE > is_focal_FALSE)
  pangolin_lineages_aggregated$parent_lineage <- unlist(lapply(
    FUN = get_parent_lineage,
    X = pangolin_lineages_aggregated$pangolin_lineage))           
  return(pangolin_lineages_aggregated)
}

# ------------------------------------------------------------------------------

#' Take sequences from the next month if insufficient samples are available from
#' current month. Assumes travel_context has complete cases with respect to iso_country
#' and date (all months have an entry for each country).
fill_travel_context_forward <- function(travel_context) {
  is_first <- T
  for (country_i in unique(travel_context$iso_country)) {
    travel_context_i <- travel_context %>% 
      filter(iso_country == country_i) %>%
      arrange(date)
    if (nrow(travel_context_i) > 1) {
      for (j in 1:(nrow(travel_context_i) - 1)) {
        n_missing <- travel_context_i[j, "n_seqs_ideal"] - travel_context_i[j, "n_seqs_accounted_for"]
        k <- j + 1
        n_extra_taken <- min(n_missing, travel_context_i[k, "n_seqs_remaining"])
        travel_context_i[k, "n_seqs_actual"] <- travel_context_i[k, "n_seqs_actual"] + n_extra_taken
        travel_context_i[j, "n_seqs_accounted_for"] <- travel_context_i[j, "n_seqs_accounted_for"] + n_extra_taken
        travel_context_i[k, "n_seqs_remaining"] <- travel_context_i[k, "n_seqs_remaining"] - n_extra_taken
        n_missing <- n_missing - n_extra_taken
      }
      if (is_first) {
        travel_context_fill_forward <- travel_context_i
        is_first <- F
      } else {
        travel_context_fill_forward <- rbind(
          travel_context_fill_forward, travel_context_i)
      }
    } else {
      travel_context_fill_forward <- travel_context
    }
  }
  return(travel_context_fill_forward)
}

#' Take sequences from the previous month if insufficient samples are available from
#' current month. Assumes travel_context has complete cases with respect to iso_country
#' and date (all months have an entry for each country).
fill_travel_context_backward <- function(travel_context) {
  is_first <- T
  for (country_i in unique(travel_context$iso_country)) {
    travel_context_i <- travel_context %>% 
      filter(iso_country == country_i) %>%
      arrange(date)
    if (nrow(travel_context_i) > 1) {
      for (j in nrow(travel_context_i):2) {
        n_missing <- travel_context_i[j, "n_seqs_ideal"] - travel_context_i[j, "n_seqs_accounted_for"]
        k <- j - 1
        n_extra_taken <- tryCatch(
          {
            min(n_missing, travel_context_i[k, "n_seqs_remaining"])
          },
          warning = function(cond) {
            message(paste(country_i, "caused a warning with n_missing:", n_missing,
                          "and n_seqs_remaining prev month:", 
                          travel_context_i[k, "n_seqs_remaining"]))
          }
        ) 
        travel_context_i[k, "n_seqs_actual"] <- travel_context_i[k, "n_seqs_actual"] + n_extra_taken
        travel_context_i[j, "n_seqs_accounted_for"] <- travel_context_i[j, "n_seqs_accounted_for"] + n_extra_taken
        travel_context_i[k, "n_seqs_remaining"] <- travel_context_i[k, "n_seqs_remaining"] - n_extra_taken
        n_missing <- n_missing - n_extra_taken
      }
      if (is_first) {
        travel_context_fill_backward <- travel_context_i
        is_first <- F
      } else {
        travel_context_fill_backward <- rbind(
          travel_context_fill_backward, travel_context_i)
      }
    } else {
      travel_context_fill_backward <- travel_context
    }
  }
  return(travel_context_fill_backward)
}

#' Report representativeness of sampling
report_travel_sampling <- function(travel_context, outdir) {
  write.table(
    x = travel_context,
    file = paste(outdir, "travel_context.txt", sep = "/"),
    row.names = F, col.names = T, quote = F, sep = "\t")
  
  travel_context_summary <- travel_context %>%
    group_by(iso_country) %>%
    summarise(actual_total_seqs = sum(n_seqs_actual),
              ideal_total_seqs = sum(n_seqs_ideal))
  
  undersampled_countries <- travel_context_summary %>% 
    filter(ideal_total_seqs > actual_total_seqs) %>%
    mutate(percent_undersampling = 1 - (actual_total_seqs / ideal_total_seqs)) %>%
    arrange(desc(percent_undersampling))
  
  if (nrow(undersampled_countries) > 0) {
    warning(
      "Some countries are underrepresented in the travel context set:\n",
      paste(capture.output(print(undersampled_countries)), collapse = "\n"))
  }
  write.table(
    x = undersampled_countries,
    file = paste(outdir, "undersampled_countries_travel_context.txt", sep = "/"),
    row.names = F, col.names = T, quote = F, sep = "\t")
}

#' Select number of travel context sequences based on estimated travel cases per 
#' month and available sequences.
get_travel_context <- function(
  n_strains, travel_cases, db_connection, outdir, qcd_gisaid_query 
) {
  # Get available sequences per country month
  context_strains <- qcd_gisaid_query %>%
    select(iso_country, date, date_str) %>%
    collect() %>%
    mutate(
      date_validity = grepl(x = date_str, pattern = "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}"),
      date = format(date, "%Y-%m-01")) %>%
    filter(date_validity) %>%
    group_by(date, iso_country) %>%
    summarize(n_seqs_available = n())
    
  travel_context <- merge(
    x = travel_cases, y = context_strains, 
    all.x = T, by = c("iso_country", "date")) %>%
    tidyr::replace_na(replace = list("n_seqs_available" = 0))
  
  # Get ideal number sequences per country month
  travel_context <- travel_context %>% mutate(
    n_seqs_ideal_raw = (n_strains * n_travel_cases) / sum(n_travel_cases),
    n_seqs_ideal = quota_largest_remainder(
      votes = travel_context$n_travel_cases, 
      n_seats = n_strains))
  print(paste("Ideally selecting", sum(travel_context$n_seqs_ideal), "travel context sequences."))
  
  # When possible, take sequences from appropriate month
  # Otherwise take from next month, and if necessary from previous month
  travel_context <- travel_context %>% mutate(
    n_seqs_actual = pmin(n_seqs_available, n_seqs_ideal),
    n_seqs_remaining = n_seqs_available - n_seqs_actual,
    n_seqs_accounted_for = n_seqs_actual)
  travel_context_fill_forward <- fill_travel_context_forward(travel_context)
  travel_context_fill_backward <- fill_travel_context_backward(travel_context_fill_forward)
  
  report_travel_sampling(travel_context_fill_backward, outdir = outdir)
  return(travel_context_fill_backward)
}

#' Sample travel context sequences.
get_travel_strains <- function(
  n_strains, travel_cases, db_connection, outdir, qcd_gisaid_query
) {
  print("Selecting number of travel context sequences based on estimated travel cases per month and available sequences")
  travel_context <- get_travel_context(
    n_strains = n_strains,
    travel_cases = travel_cases,
    db_connection = db_connection,
    outdir = outdir,
    qcd_gisaid_query = qcd_gisaid_query)
  if (sum(travel_context$n_seqs_actual) == 0) {
    warning("No travel context sequences selected. If this is not desired, either increase timeframe or travel context scale factor.")
    travel_strains <- data.frame(
      gisaid_epi_isl = NA,
      iso_country = NA,
      date = NA,
      date_str = NA,
      pangolin_lineage = NA)
    write.csv(x = travel_strains,
            file = paste(outdir, "travel_strains.txt", sep = "/"),
            row.names = F)
    return(travel_strains)
  }
  
  print("Getting available sequences per country month from database table gisaid_sequence.")
  context_strains <- qcd_gisaid_query %>%
    select(gisaid_epi_isl, iso_country, date, date_str, pangolin_lineage) %>%
    collect() %>%
    mutate(date_validity = grepl(x = date_str, pattern = "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}"), 
           date = format(date, "%Y-%m-01")) %>%
    filter(date_validity, iso_country %in% travel_context$iso_country)
  
  travel_context_nonzero_samples <- travel_context %>% filter(n_seqs_actual > 0)
  n_country_months <- nrow(travel_context_nonzero_samples)
  l <- 0
  for (iso_country_i in unique(travel_context_nonzero_samples$iso_country)) {
    travel_context_i <- travel_context_nonzero_samples %>% filter(iso_country == iso_country_i)
    context_strains_i <- context_strains %>% filter(iso_country == iso_country_i)
    for (j in 1:nrow(travel_context_i)) {
      date_k <- travel_context_i[j, "date"]
      n_seqs_actual_k <- travel_context_i[j, "n_seqs_actual"]
      status <- paste("sampling", n_seqs_actual_k, "travel context strains from", iso_country_i, "in", date_k)
      l <- l + 1
      if (l %% 10 == 0) {
        print(paste(round(l * 100 / n_country_months), 
                    "% done... currently", status))
      }
      travel_strains_k <- tryCatch(
        {
          context_strains_i %>%
            filter(date == date_k) %>%
            sample_n(size = n_seqs_actual_k, replace = F)
        },
        error = function(cond) {
          message(paste(status, "resulted in an error with message", cond))
          stop()
        }
      )
      if (l == 1) {
        travel_strains <- travel_strains_k
        is_first <- F
      } else {
        travel_strains <- rbind(travel_strains, travel_strains_k)
      }
    }
  }
  
  write.csv(x = travel_strains,
            file = paste(outdir, "travel_strains.txt", sep = "/"),
            row.names = F)
  return(travel_strains)
}

# ------------------------------------------------------------------------------

#' Run nextstrain priority script to rank context sequences by genetic proximity
#' to focal sequenes. Writes out file with strains and corresponding priority.
#' @return table with strains and corresponding priority
run_nextstrain_priority <- function(
  focal_strains, nonfocal_strains, prefix, outdir, 
  python_path, priorities_script, reference, verbose = F
) {
  n_focal_strains <- length(focal_strains)
  n_nonfocal_strains <- length(nonfocal_strains)
  
  if (n_focal_strains == 0) {
    warning(paste("No focal sequences found for", prefix))
    return(NA)
  } else if (n_nonfocal_strains == 0) {
    warning(paste("No context sequences found for", prefix))
    return(NA)
  } else if (n_nonfocal_strains > 400000) {
    stop(paste("Too many sequences for priorities.py:", 
               n_nonfocal_strains, "non-focal sequences."))
  }
  print(paste("Running priorities.py for", prefix, "with", 
              length(focal_strains), "focal sequences and", 
              length(nonfocal_strains), "non-focal sequences."))

  outfile_focal_strains <- paste(
    outdir, "/", prefix, "_focal_strains.txt", sep = "")
  outfile_context_strains <- paste(
    outdir, "/", prefix, "_nonfocal_strains.txt", sep = "")
  write.table(
    x = c("strain", focal_strains), row.names = F, quote = F, col.names = F,
    file = outfile_focal_strains)
  write.table(
    x = c("strain", nonfocal_strains), row.names = F, quote = F, 
    col.names = F, file = outfile_context_strains)
  
  outfile_priorities <- paste(outdir, "/", prefix, "_priorities.txt", sep = "")
  priorities_command <- paste(
    python_path, priorities_script,
    "--focal-strains", outfile_focal_strains,
    "--context-strains", outfile_context_strains,
    "--reference", reference,
    "--outfile", outfile_priorities,
    "--config-filepath", "workdir/input/config.yml",
    "--automated")
  if (verbose) {
    print(priorities_command)
  }
  system(command = priorities_command)
  
  priorities <- read.delim(
    file = paste(outdir, "/", prefix, "_priorities.txt", sep = ""),
    header = F,
    col.names = c("strain", "priority", "focal_strain"),
    stringsAsFactors = F) %>%
    arrange(desc(priority))
  
  system(command = paste("rm", outfile_focal_strains))
  system(command = paste("rm", outfile_context_strains))
  # system(command = paste("rm", outfile_priorities))
  
  return(priorities)
}

get_similarity_strains <- function(
  db_connection, qcd_gisaid_query, similarity_context_scale_factor, lineages,
  priorities_script = "database/python/priorities_from_database.py",
  python_path, reference, verbose = F, unique_context_only, focal_country
) {
  n_lineages <- nrow(lineages)
  initiated_df <- F
  for (i in 1:n_lineages) {
    lineage <- lineages$pangolin_lineage[i]
    lineages_included_str <- lineages$lineages_aggregated[i]
    lineages_included <- strsplit(lineages_included_str, split = ", ")[[1]]
    
    focal_seqs <- qcd_gisaid_query %>%
      filter(country == focal_country,
             pangolin_lineage %in% !! lineages_included) %>%
      select(strain) %>%
      collect()
    prospective_context_seqs <- get_prospective_context_seqs(
      qcd_gisaid_query = qcd_gisaid_query,
      db_connection = db_connection,
      unique_context_only = unique_context_only,
      lineages_included = lineages_included,
      focal_country = focal_country)
    
    if (nrow(focal_seqs) == 0 | nrow(prospective_context_seqs) == 0) {
      print(paste("No focal sequences from lineage", lineage))
      next
    }
    
    focal_strains <- focal_seqs$strain
    prospective_context_strains <- prospective_context_seqs$strain
    priorities <- run_nextstrain_priority(
      focal_strains = focal_strains, 
      nonfocal_strains = prospective_context_strains,
      prefix = lineage,
      outdir = outdir, 
      priorities_script = priorities_script,
      python_path = python_path, 
      reference = reference, 
      verbose = verbose)
    
    n_similarity_seqs <- ceiling(length(focal_strains) * similarity_context_scale_factor)
    similarity_strains_i <- priorities %>%
        arrange(desc(priority)) %>%
        mutate(priority_idx = 1:n()) %>%
        top_n(n = min(n_similarity_seqs, n()), wt = -priority_idx) %>%
        select(strain) %>%
        mutate(pangolin_lineage = lineage) %>%
        left_join(prospective_context_seqs, by = "strain")
    
    if (!initiated_df) {
      similarity_strains <- similarity_strains_i
      initiated_df <- T
    } else {
      similarity_strains <- rbind(similarity_strains, similarity_strains_i)
    }
  }
  return(similarity_strains)
}

#' Filter context set to the earliest sequences with unique AA mutations.
get_prospective_context_seqs <- function(
  qcd_gisaid_query, db_connection, unique_context_only, lineages_included, focal_country
) {
  prospective_context_seqs <- qcd_gisaid_query %>%
      filter(country != focal_country,
             pangolin_lineage %in% !! lineages_included) %>%
      select(strain, gisaid_epi_isl, country, date) %>%
      collect()

  if (unique_context_only) {
    print("Retaining only earliest of prospective context strains with identical amino acid mutations.")
    # Generate string of prospective context strains suitable for query
    quoted_strains = lapply(prospective_context_seqs$gisaid_epi_isl, function(elt) {
        DBI::dbQuoteString(conn = db_connection, elt)
    })
    strains_sql = paste0('(', do.call(paste, c(quoted_strains, sep = ',')), ')')
    
    # Select earliest of strains with unique amino acid mutations
    prospective_context_seqs_sql <- paste("select",
        "min(strain order by date) as strain,",
        "min(gisaid_epi_isl order by date) as gisaid_epi_isl,",
        "min(country order by date) as country,",
        "min(date) as date,",
        "count(*) as n_other_seqs_w_same_aa_mutations",
      "from (",
        "select",
          "gs.strain,",
          "gs.gisaid_epi_isl,",
          "gs.country,",
          "string_agg(aa_mutation, ',' order by aa_mutation) as aa_mutations,",
          "gs.date",
        "from",
          "gisaid_api_sequence gs",
          "left join gisaid_api_sequence_nextclade_mutation_aa gsnma on gs.gisaid_epi_isl = gsnma.gisaid_epi_isl",
        "where",
          "gs.gisaid_epi_isl in", strains_sql,
        "group by gs.gisaid_epi_isl ) seqs_w_aa_mutations",
      "group by aa_mutations;")
    
   prospective_context_seqs <- DBI::dbGetQuery(
      conn = db_connection, statement = prospective_context_seqs_sql)
  }
  return(prospective_context_seqs)
}

# ------------------------------------------------------------------------------

#' Write out an alignment for each lineage containing all QC-passed focal seqs,
#' travel & genetic similarity context seqs, and outgroup seqs.
#' @return Dataframe giving size of each alignment.
write_out_alignments <- function(
  lineages, travel_strains, similarity_strains, outgroup_gisaid_epi_isls, outdir,
  qcd_gisaid_query, db_connection, mask_from_start, mask_from_end, focal_country
) {
  alignments_generated <- c()
  n_lineages <- nrow(lineages)
  initiated_df <- F
  for (i in 1:n_lineages) {
    lineage <- lineages$pangolin_lineage[i]
    lineages_included_str <- lineages$lineages_aggregated[i]
    lineages_included <- strsplit(lineages_included_str, split = ", ")[[1]]
    
    focal_strains_i <- qcd_gisaid_query %>%
      filter(country == focal_country,
             pangolin_lineage %in% !! lineages_included) %>%
      collect()
    
    if (nrow(focal_strains_i) == 0) {
      print(paste("No focal strains for lineage", lineage, "so not writing out an alignment."))
      next
    }
    
    travel_strains_i <- travel_strains %>% 
      filter(pangolin_lineage %in% lineages_included)
    similarity_strains_i <- similarity_strains %>% 
      filter(pangolin_lineage %in% lineages_included)
    alignment_strains_i <- unique(c(
      outgroup_gisaid_epi_isls, 
      focal_strains_i$gisaid_epi_isl,
      travel_strains_i$gisaid_epi_isl, 
      similarity_strains_i$gisaid_epi_isl))
    
    cat(
      "Lineage", lineage, "alignment has", length(alignment_strains_i), "sequences:\n", 
      nrow(focal_strains_i), "focal seqs\n",
      nrow(travel_strains_i), "travel context seqs\n",
      nrow(similarity_strains_i), "genetic similarity context seqs and\n",
      length(outgroup_gisaid_epi_isls), "outgroup sequences\n")
    
    metadata_i <- dplyr::tbl(db_connection, "gisaid_api_sequence") %>%
      filter(gisaid_epi_isl %in% !! alignment_strains_i) %>%
      select(gisaid_epi_isl, strain, date, region_original, country, division,
             pangolin_lineage, originating_lab, submitting_lab, authors) %>%
      collect() %>%
      mutate(
        tree_pangolin_lineage = case_when(
          gisaid_epi_isl %in% outgroup_gisaid_epi_isls ~ pangolin_lineage,
          T ~ lineage),
        tree_label = paste(gisaid_epi_isl, date, sep = "|"),
        travel_context = case_when(
          gisaid_epi_isl %in% travel_strains_i$gisaid_epi_isl ~ T,
          T ~ F),
        similarity_context = case_when(
          gisaid_epi_isl %in% similarity_strains_i$gisaid_epi_isl ~ T,
          T ~ F),
        focal_sequence = case_when(
          gisaid_epi_isl %in% focal_strains_i$gisaid_epi_isl ~ T,
          T ~ F))
    missing_strains <- alignment_strains_i[!(alignment_strains_i %in% metadata_i$gisaid_epi_isl)]
    if (length(missing_strains) > 0) {
      stop(paste(paste0(missing_strains, collapse = ", "), "selected strains not found in table gisaid_sequence."))
    }
    col_order <- c("tree_label", colnames(metadata_i)[colnames(metadata_i) != "tree_label"])  # reorder for figtree
    write.table(
      x = metadata_i %>% select(all_of(col_order)),
      file = paste(outdir, "/", lineage, "_metadata.tsv", sep = ""),
      sep = "\t",
      row.names = F)

    header_mapping <- metadata_i$tree_label
    names(header_mapping) <- metadata_i$gisaid_epi_isl
      
    export_seqs_as_fasta(
      db_connection = db_connection, 
      sample_names = alignment_strains_i,
      seq_outfile = paste(outdir, "/", lineage, ".fasta", sep = ""),
      table = "gisaid_api_sequence",
      sample_name_col = "gisaid_epi_isl",
      seq_col = "seq_aligned",
      header_mapping = header_mapping,
      mask_from_start = mask_from_start,
      mask_from_end = mask_from_end)
    
    alignment_info_i <- data.frame(
      lineage = lineage, n_seqs = length(alignment_strains_i))
    if (!initiated_df) {
      alignment_info <- alignment_info_i
      initiated_df <- T
    } else {
      alignment_info <- rbind(alignment_info, alignment_info_i)
    }
  }
  return(alignment_info %>% arrange(desc(n_seqs)))
}
  