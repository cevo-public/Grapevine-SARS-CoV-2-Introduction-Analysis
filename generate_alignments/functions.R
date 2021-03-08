#' @param qcd_gisaid_query Query of database table gisaid_sequence with QC filters
#' @return qcd_gisaid_query Query expanded to include only selected swiss samples
downsample_swiss_sequences <- function (
  qcd_gisaid_query, db_connection, max_sampling_frac
) {
  sampling_data <- get_weekly_case_and_seq_data(
    db_connection = db_connection, qcd_gisaid_query = qcd_gisaid_query) %>%
    mutate(
      n_seqs_viollier = tidyr::replace_na(n_seqs_viollier, replace = 0),
      n_seqs_other = tidyr::replace_na(n_seqs_other, replace = 0),
      n_seqs = n_seqs_viollier + n_seqs_other,
      max_seqs = floor(n_conf_cases * max_sampling_frac),
      n_to_sample = pmin(n_seqs, max_seqs))
  
  qcd_gisaid_query_temp <- qcd_gisaid_query %>%
    filter(country == "Switzerland") %>%
    mutate(week = date_trunc('week', date))
  
  sampled_strains <- c()
  for (i in 1:nrow(sampling_data)) {
    week_i <- sampling_data[[i, "week"]]
    n_samples_i <- as.numeric(sampling_data[[i, "n_to_sample"]])
    n_actual_samples <- nrow(qcd_gisaid_query_temp %>%
                               filter(week == week_i) %>%
                               select(strain) %>%
                               collect())
    print(paste("sampling", n_samples_i, "Swiss samples from", n_actual_samples, "in week", week_i))
    
    sampled_strains_i <- unname(unlist(qcd_gisaid_query_temp %>%
      filter(week == week_i) %>%
      select(strain) %>%
      collect() %>%
      sample_n(size = n_samples_i, replace = F)))
    
    sampled_strains <- c(sampled_strains, sampled_strains_i)
  }

  qcd_gisaid_query_swiss_downsampled <- qcd_gisaid_query %>%
    filter(!(country == 'Switzerland' & !(strain %in% !! sampled_strains)))
  
  return(qcd_gisaid_query_swiss_downsampled)
}

# ------------------------------------------------------------------------------

#' Get number of cases by exposure country and month of case confirmation from 
#' BAG meldeformular.
#' Null entries for 'exp_land' or entries of 'Schweiz' in BAG meldeformular are not reported.
get_exposures_per_country_month <- function(db_connection, min_date, max_date) {
  print("Getting number of cases by exposure country and month of case confirmation from BAG meldeformular.")
  exposures_per_country_month <- dplyr::tbl(
    db_connection, "bag_meldeformular") %>%
    filter(!is.na(exp_land), exp_land != "Schweiz", 
           fall_dt <= !! max_date, fall_dt >= !! min_date) %>%
    select(exp_land, fall_dt) %>%
    collect() %>%
    mutate(date = format(fall_dt, "%Y-%m-01")) %>%
    group_by(date, exp_land) %>%
    summarise(n_exposures = n()) %>%
    mutate(exp_land = recode(
      exp_land, 
      "Deutschland ohne nähere Angaben" = "Deutschland"))
  country_translation <- dplyr::tbl(db_connection, "country") %>% 
    collect() %>%
    rename("exp_land" = "german_name")
  exposures_per_country_month_2 <- left_join(
    x = exposures_per_country_month, y = country_translation, 
    na_matches = "never", by = "exp_land") %>%
    select(-c(english_name)) %>%
    rename("iso_code" = "iso3166_alpha3_code")
  return(exposures_per_country_month_2)
}

#' Estimate number of tourists arriving in Switzerland by origin country and 
#' month of arrival. 
#' @param db_connection
#' @param max_date Character date, e.g. "2021-01-30" to extrapolate data to.
get_tourist_arrivals_per_country_month <- function(db_connection, min_date, max_date) {
  print("Estimating number of tourists arriving in Switzerland by origin country and month of arrival.")
  tourist_arrivals_per_country_month <- dplyr::tbl(
    db_connection, "ext_fso_tourist_accommodation") %>%
    filter(date <= !! max_date, date >= !! min_date) %>%
    collect() %>%
    tidyr::complete(
      date = seq.Date(from = min(date), to = as.Date(max_date), by = "month"), 
      iso_code) %>%  
    group_by(iso_code) %>%
    arrange(date) %>%
    tidyr::fill(n_arrivals) %>% # fill in missing months with prev month's value
    rename(n_tourist_arrivals = n_arrivals) %>%
    filter(!is.na(iso_code))
  return(tourist_arrivals_per_country_month)
}

# Get number of cross-border commuter permits for Switzerland by origin country 
# and month.
#' @param db_connection
#' @param max_date Character date, e.g. "2021-01-30" to extrapolate data to.
get_commuter_permits_per_country_month <- function(db_connection, min_date, max_date) { 
  print("Getting number of cross-border commuter permits for Switzerland by origin country and month.")
  commuter_permits_per_country_month <- dplyr::tbl(
    db_connection, "ext_fso_cross_border_commuters") %>%
    filter(date <= !! max_date, date >= !! min_date) %>%
    collect() %>%
    tidyr::complete(
      date = seq.Date(from = as.Date(min_date), to = as.Date(max_date), by = "month"), 
      iso_code, wirtschaftsabteilung) %>% 
    group_by(iso_code, wirtschaftsabteilung) %>%
    arrange(date) %>%
    tidyr::fill(n_permits) %>% # fill in months with quarterly values, missing quarters with previous quarter's value
    group_by(iso_code, date) %>%
    summarise(n_commuter_permits = sum(n_permits)) %>%  # sum across sectors
    filter(!is.na(iso_code))
  return(commuter_permits_per_country_month)
}

#' Get average daily infectious population (per million) by country and month.
#' @param db_connection
#' @param max_date Character date, e.g. "2021-01-30" to extrapolate data to.
get_avg_infectious_per_country_month <- function(db_connection, min_date, max_date) {
  print("Geting average daily infectious population (per million) by country and month.")
  avg_infectious_per_country_month <- dplyr::tbl(
    db_connection, "ext_owid_global_cases") %>%
    filter(date <= !! max_date, date >= !! min_date) %>%
    select(iso_code, date, new_cases_per_million) %>%
    group_by(iso_code) %>%
    arrange(date) %>%
    mutate(cumul_cases_per_million = cumsum(new_cases_per_million)) %>%
    collect() %>%
    mutate(n_infectious_cases_per_million = lead(cumul_cases_per_million, n = 10) - cumul_cases_per_million) %>%
    mutate(date = as.Date(format(date, "%Y-%m-01"))) %>%
    group_by(date, iso_code) %>%
    summarise(avg_daily_n_infectious_per_million = mean(n_infectious_cases_per_million, na.rm = T)) %>%
    tidyr::complete(
      date = seq.Date(from = min(date), to = as.Date(max_date), by = "month"), 
      avg_daily_n_infectious_per_million) %>% 
    arrange(date) %>%
    tidyr::fill(avg_daily_n_infectious_per_million) %>% # fill in months with previous monthly value
    filter(!is.na(iso_code))
  return(avg_infectious_per_country_month)  
}

#' Estimate number of infectious arrivals into Switzerland by origin country and
#' month of arrival.
#' @param db_connection
#' @param max_date Character date, e.g. "2021-01-30" to extrapolate data to.
get_infected_arrivals_per_country_month <- function(db_connection, min_date, max_date, outdir) {
  print("Estimating number of infectious arrivals into Switzerland by origin country and month of arrival.")
  tourists <- get_tourist_arrivals_per_country_month(db_connection, min_date, max_date)
  commuters <- get_commuter_permits_per_country_month(db_connection, min_date, max_date)
  infectious <- get_avg_infectious_per_country_month(db_connection, min_date, max_date)
  arrivals <- merge(x = tourists, y = commuters, all = T, 
                    by = c("date", "iso_code")) %>%
    tidyr::replace_na(replace = list("n_tourist_arrivals" = 0, 
                                     "n_commuter_permits" = 0)) %>%
    mutate(n_arrivals = n_tourist_arrivals + n_commuter_permits)
  infectious_arrivals <- merge(x = arrivals, y = infectious, 
                               by = c("date", "iso_code"), all.x = T)
  warn_missing_infectious_arrivals(infectious_arrivals)
  infectious_arrivals <- infectious_arrivals %>% 
    filter(!is.na(avg_daily_n_infectious_per_million)) %>%
    mutate(n_infectious_arrivals = avg_daily_n_infectious_per_million * n_arrivals / 1E6)
  report_infectious_arrivals(infectious_arrivals, outdir)
  return(infectious_arrivals)
}

warn_missing_infectious_arrivals <- function(infectious_arrivals) {
  missing_case_data <- infectious_arrivals %>%
    filter(is.na(avg_daily_n_infectious_per_million), 
           date >= as.Date("2020-03-01")) %>%
    group_by(iso_code, date) %>%
    summarize(missing_case_data = T) 
  if (nrow(missing_case_data) > 0) {
    warning(
      "These countries are missing case data from after 2020-03-01 but have tourist/commuter arrivals into Switzerland. These arrivals will be ignored for context sequence subsampling.",
      paste(capture.output(print(missing_case_data)), collapse = "\n"))
  }
}

report_infectious_arrivals <- function(infectious_arrivals, outdir) {
  print("Plotting estimated infectious arrivals per country and month.")
  require(ggplot2)
  write.table(
    x = infectious_arrivals,
    file = paste(outdir, "estimated_monthly_infectious_arrivals.txt", sep = "/"),
    row.names = F, col.names = T, quote = F, sep = "\t")
  top_10_countries <- infectious_arrivals %>%
    group_by(iso_code) %>%
    summarize(total_infectious_arrivals = sum(n_infectious_arrivals)) %>%
    top_n(n = 8, wt = total_infectious_arrivals) %>%
    select(iso_code)
  p <- ggplot(
    data = infectious_arrivals %>% 
      mutate(source_country = case_when(
        iso_code %in% top_10_countries$iso_code ~ iso_code,
        T ~ "Other")),
    aes(x = format(date, "%Y-%m"), y = n_infectious_arrivals)) +
    geom_col(aes(fill = source_country)) + 
    theme_bw() + 
    labs(x = element_blank(), y = "No. individuals") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  ggsave(
    filename = paste(outdir, "estimated_monthly_infectious_arrivals.png", sep = "/"), 
    plot = p)
}

#' Get a set of strains proportional to estimated travel connections to Switzerland.
#' Estimated infectious individuals arriving and travel exposure cases are both
#' in absolute numbers, so sum for total number infectious inidiviuals arriving
#' in Switzerland from each country.
#' @param db_connection
get_travel_cases <- function(db_connection, min_date, max_date, outdir) {
  print("Pulling tourist accommodation and cross-border commuter stats from FSO into the database.")
  source("database/R/import_fso_travel_stats.R")
  import_fso_tourist_accommodation(db_connection = db_connection)
  import_fso_cross_border_commuters(db_connection = db_connection)
  
  print("Querying tourist and exposure information from database.")
  infected_arrivals <- get_infected_arrivals_per_country_month(db_connection, min_date, max_date, outdir)
  exposures <- get_exposures_per_country_month(db_connection, min_date, max_date)
  
  print("Merging arrival and exposure data.")
  travel_cases <- merge(
    x = infected_arrivals, y = exposures, 
    all = T, by = c("date", "iso_code")) %>%
    tidyr::replace_na(
      replace = list("n_exposures" = 0, "n_infectious_arrivals" = 0)) %>%
    mutate(n_travel_cases = n_exposures + n_infectious_arrivals) %>%
    arrange(date, desc(n_infectious_arrivals))
  
  to_remove <- travel_cases %>%
    filter(is.na(iso_code)) %>%
    select(origin, exp_land, n_arrivals, n_exposures)
  if (length(to_remove) > 0) {
    warning("Not adding travel context sequences for these entries because no valid iso country code found.\n",
            paste(capture.output(print(to_remove)), collapse = "\n"))
    travel_cases <- travel_cases %>% filter(!is.na(iso_code))
  }
  
  print("Completing travel case data with zero values for months with no travel cases.")
  travel_cases <- travel_cases %>%
    tidyr::complete(
      date = seq.Date(from = min(date), to = as.Date(max_date), by = "month"), 
      iso_code) %>%  
    group_by(iso_code) %>%
    arrange(date) %>%
    tidyr::replace_na(replace = list("n_travel_cases" = 0))  # fill in missing months with 0 value
  
  print(paste("Writing out results to", outdir))
  write.csv(
    x = travel_cases, 
    file = paste(outdir, "estimated_travel_cases.csv", sep = "/"), 
    row.names = F)
  
  return(travel_cases)
}

# ------------------------------------------------------------------------------

#' Get the parent lineage for a pangolin lineage, or return NA if there is no 
#' valid parent lineage.
get_parent_lineage <- function(pangolin_lineage) {
  if (grepl(x = pangolin_lineage, pattern = "\\.")) {
    # Given B.1.1.1, return B.1.1
    return(sub(".[^.]+$", "", pangolin_lineage))
  } else if (nchar(pangolin_lineage) == 1 && utf8ToInt(pangolin_lineage) %in% 66:90) {  
    # According to pangolin website, C is an alias for B.1.1.1
    # Given C, return B.1.1.1
    lineage_letter <- intToUtf8(utf8ToInt(pangolin_lineage) - 1)
    return(paste(lineage_letter, ".1.1.1", sep = ""))
  } else {
    warning(paste("Cannot find a valid parent lineage for", pangolin_lineage, "\n"))
    return(NA)
  }
}

#' For pangolin lineages where > 50% of samples are Swiss, assign samples to the 
#' parent lineage. The goal here is that the MRCA of all lineages is outside of 
#' Switzerland. 
aggregate_predominatly_swiss_lineages <- function(pangolin_lineages) {
  pangolin_lineages_aggregated <- pangolin_lineages %>% 
    mutate(
      pangolin_lineage_rollup = case_when(
        is_swiss_TRUE > is_swiss_FALSE & !is.na(parent_lineage) ~ parent_lineage,
        T ~ pangolin_lineage),
      is_high_level_swiss_lineage = should_aggregate & is.na(parent_lineage)) %>%
    group_by(pangolin_lineage_rollup) %>%
    summarize(is_swiss_TRUE = sum(is_swiss_TRUE),
              is_swiss_FALSE = sum(is_swiss_FALSE),
              n_lineages_aggregated = sum(n_lineages_aggregated),
              lineages_aggregated = paste0(lineages_aggregated, collapse = ", "))
  
  pangolin_lineages_aggregated <- pangolin_lineages_aggregated %>%
    mutate(pangolin_lineage = pangolin_lineage_rollup,
           should_aggregate = is_swiss_TRUE > is_swiss_FALSE)
  pangolin_lineages_aggregated$parent_lineage <- unlist(lapply(
    FUN = get_parent_lineage,
    X = pangolin_lineages_aggregated$pangolin_lineage
  ))           
  return(pangolin_lineages_aggregated)
}

#' QC gisaid data and split by pangolin lineage, aggregating lineages that are 
#' predominantly Swiss into the parent lineage.
get_pangolin_lineages <- function(db_connection, outdir, qcd_gisaid_query) {
  print("Querying database for pangolin lineages.")
  pangolin_lineages <- qcd_gisaid_query %>%
    mutate(is_swiss = country == "Switzerland") %>%
    group_by(pangolin_lineage, is_swiss) %>%
    summarize(n_seqs = n()) %>%
    collect() %>%
    tidyr::pivot_wider(
      names_from = is_swiss, names_prefix = "is_swiss_", values_from = n_seqs) %>%
    tidyr::replace_na(replace = list("is_swiss_TRUE" = 0, "is_swiss_FALSE" = 0)) %>%  # because by default 0s show up as NA after group_by %>% summarize
    filter(is_swiss_TRUE > 0) %>%
    mutate(
      should_aggregate = is_swiss_TRUE > is_swiss_FALSE,
      parent_lineage = get_parent_lineage(pangolin_lineage),
      n_lineages_aggregated = 1,
      lineages_aggregated = pangolin_lineage)
  print("Aggregating predominantly Swiss lineages into parent lineage.")
  pangolin_lineages_aggregated <- pangolin_lineages
  while(any(!is.na(pangolin_lineages_aggregated$parent_lineage) & 
            pangolin_lineages_aggregated$should_aggregate)) {
    pangolin_lineages_aggregated <- aggregate_predominatly_swiss_lineages(
      pangolin_lineages_aggregated)
  }
  write.csv(
    x = pangolin_lineages_aggregated,
    file = paste(outdir, "pangolin_lineages_aggregated.csv", sep = "/"),
    row.names = F)
  return(pangolin_lineages_aggregated)
}

# ------------------------------------------------------------------------------

#' Take sequences from the next month if insufficient samples are available from
#' current month. Assumes travel_context has complete cases with respect to iso_code
#' and date (all months have an entry for each country).
fill_travel_context_forward <- function(travel_context) {
  is_first <- T
  for (country_i in unique(travel_context$iso_code)) {
    travel_context_i <- travel_context %>% 
      filter(iso_code == country_i) %>%
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
#' current month. Assumes travel_context has complete cases with respect to iso_code
#' and date (all months have an entry for each country).
fill_travel_context_backward <- function(travel_context) {
  is_first <- T
  for (country_i in unique(travel_context$iso_code)) {
    travel_context_i <- travel_context %>% 
      filter(iso_code == country_i) %>%
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
    group_by(iso_code) %>%
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
    select(country, date, date_str) %>%
    collect() %>%
    mutate(
      date_validity = grepl(x = date_str, pattern = "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}"),
      date = format(date, "%Y-%m-01"),
           iso_code = countrycode::countrycode(
             sourcevar = country, origin = "country.name", 
             destination = "iso3c")) %>%
    filter(date_validity) %>%
    group_by(date, country, iso_code) %>%
    summarize(n_seqs_available = n()) %>% 
    filter(!is.na(iso_code))
    
  travel_context <- merge(
    x = travel_cases, y = context_strains, 
    all.x = T, by = c("iso_code", "date")) %>%
    tidyr::replace_na(replace = list("n_seqs_available" = 0))
  
  # Get ideal number sequences per country month
  travel_context <- travel_context %>% mutate(
    n_seqs_ideal = round((n_strains * n_travel_cases) / sum(n_travel_cases)))
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
  
  print("Getting available sequences per country month from database table gisaid_sequence.")
  context_strains <- qcd_gisaid_query %>%
    select(gisaid_epi_isl, country, date, date_str, pangolin_lineage) %>%
    collect() %>%
    mutate(date_validity = grepl(x = date_str, pattern = "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}"), 
           date = format(date, "%Y-%m-01"),
           iso_code = countrycode::countrycode(
             sourcevar = country, origin = "country.name", 
             destination = "iso3c")) %>%
    filter(date_validity, iso_code %in% travel_context$iso_code)
  
  travel_context_nonzero_samples <- travel_context %>% filter(n_seqs_actual > 0)
  n_country_months <- nrow(travel_context_nonzero_samples)
  l <- 0
  for (iso_code_i in unique(travel_context_nonzero_samples$iso_code)) {
    travel_context_i <- travel_context_nonzero_samples %>% filter(iso_code == iso_code_i)
    context_strains_i <- context_strains %>% filter(iso_code == iso_code_i)
    for (j in 1:nrow(travel_context_i)) {
      date_k <- travel_context_i[j, "date"]
      n_seqs_actual_k <- travel_context_i[j, "n_seqs_actual"]
      status <- paste("sampling", n_seqs_actual_k, "travel context strains from", iso_code_i, "in", date_k)
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
  python_path, priorities_script, reference
) {
  n_focal_strains <- length(focal_strains)
  n_nonfocal_strains <- length(nonfocal_strains)
  
  if (n_focal_strains == 0) {
    warning(paste("No focal sequences found for", prefix))
    return(NA)
  } else if (n_nonfocal_strains == 0) {
    warning(paste("No context sequences found for", prefix))
    return(NA)
  } else if (n_nonfocal_strains > 1000000) {
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
    "--automated")
  system(command = priorities_command)
  
  priorities <- read.delim(
    file = paste(outdir, "/", prefix, "_priorities.txt", sep = ""),
    header = F,
    col.names = c("strain", "priority"),
    stringsAsFactors = F) %>%
    arrange(desc(priority))
  
  system(command = paste("rm", outfile_focal_strains))
  system(command = paste("rm", outfile_context_strains))
  system(command = paste("rm", outfile_priorities))
  
  return(priorities)
}

get_similarity_strains <- function(
  db_connection, qcd_gisaid_query, similarity_context_scale_factor, lineages,
  priorities_script = "database/python/priorities_from_database.py",
  python_path, reference
) {
  n_lineages <- nrow(lineages)
  for (i in 1:n_lineages) {
    lineage <- lineages$pangolin_lineage[i]
    lineages_included_str <- lineages$lineages_aggregated[i]
    lineages_included <- strsplit(lineages_included_str, split = ", ")[[1]]
    
    focal_seqs <- qcd_gisaid_query %>%
      filter(country == "Switzerland", 
             pangolin_lineage %in% !! lineages_included) %>%
      select(strain) %>%
      collect()
    prospective_context_seqs <- qcd_gisaid_query %>%
      filter(country != "Switzerland",
             pangolin_lineage %in% !! lineages_included) %>%
      select(strain, gisaid_epi_isl, country, date_str, date) %>%
      collect()
    
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
      reference = reference)
    
    n_similarity_seqs <- ceiling(length(focal_strains) * similarity_context_scale_factor)
    similarity_strains_i <- priorities %>%
        arrange(desc(priority)) %>%
        mutate(priority_idx = 1:n()) %>%
        top_n(n = min(n_similarity_seqs, n()), wt = priority_idx) %>%
        select(strain) %>%
        mutate(pangolin_lineage = lineage) %>%
        left_join(prospective_context_seqs, by = "strain")
    
    if (i == 1) {
      similarity_strains <- similarity_strains_i
    } else {
      similarity_strains <- rbind(similarity_strains, similarity_strains_i)
    }
  }
  return(similarity_strains)
}

# ------------------------------------------------------------------------------

#' Write out an alignment for each lineage containing all QC-passed Swiss seqs,
#' travel & genetic similarity context seqs, and outgroup seqs.
#' @return Dataframe giving size of each alignment.
write_out_alignments <- function(
  lineages, travel_strains, similarity_strains, outgroup_gisaid_epi_isls, outdir,
  qcd_gisaid_query, db_connection
) {
  alignments_generated <- c()
  n_lineages <- nrow(lineages)
  for (i in 1:n_lineages) {
    lineage <- lineages$pangolin_lineage[i]
    lineages_included_str <- lineages$lineages_aggregated[i]
    lineages_included <- strsplit(lineages_included_str, split = ", ")[[1]]
    
    focal_strains_i <- qcd_gisaid_query %>%
      filter(country == "Switzerland", 
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
    
    metadata_i <- dplyr::tbl(db_connection, "gisaid_sequence") %>%
      filter(gisaid_epi_isl %in% !! alignment_strains_i) %>%
      select(gisaid_epi_isl, strain, date_str, date, region, country, division, 
             country_exposure, nextstrain_clade, pangolin_lineage, 
             originating_lab, submitting_lab, authors) %>%
      collect() %>%
      mutate(
        tree_pangolin_lineage = case_when(
          gisaid_epi_isl %in% outgroup_gisaid_epi_isls ~ pangolin_lineage,
          T ~ lineage),
        tree_label = paste(gisaid_epi_isl, date_str, sep = "|"),
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
    write.csv(
      x = metadata_i, 
      file = paste(outdir, "/", lineage, "_metadata.csv", sep = ""),
      row.names = F)
    
    header_mapping <- metadata_i$tree_label
    names(header_mapping) <- metadata_i$gisaid_epi_isl
      
    export_seqs_as_fasta(
      db_connection = db_connection, 
      sample_names = alignment_strains_i,
      seq_outfile = paste(outdir, "/", lineage, ".fasta", sep = ""),
      table = "gisaid_sequence",
      sample_name_col = "gisaid_epi_isl",
      seq_col = "aligned_seq",
      header_mapping = header_mapping)
    
    alignment_info_i <- data.frame(
      lineage = lineage, n_seqs = length(alignment_strains_i))
    if (i == 1) {
      alignment_info <- alignment_info_i
    } else {
      alignment_info <- rbind(alignment_info, alignment_info_i)
    }
  }
  return(alignment_info %>% arrange(desc(n_seqs)))
}
  