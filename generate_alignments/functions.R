#' Downsample Swiss sequences on GISAID proportionally to weekly confirmed cases.
#' Optionally downsample proportionally to weekly confirmed cases in each canton.
#' If there aren't enough sequences in a week from the region, all sequences are taken.
#' When confirmed cases are not attributed to a region, the proportional number of sequences are taken without regard to region.
#' Optionally favor samples with recorded foreign exposure information per BAG meldeformular or GISAID metadata before filling up the week's sequence quota with non-exposed/unknown samples.
#' @param qcd_gisaid_query Query of database table gisaid_sequence with QC filters.
#' @param max_sampling_frac Take no more than this fraction of confirmed cases each week.
#' @param favor_exposures If true, first take sequences with recorded foreign exposures.
#' @param subsample_by_canton If true, take sequences proportional to confirmed cases per canton.
#' @param output_summary_table If true, output table summarizing sampled sequences per region per week compare to confirmed cases 
#' @return qcd_gisaid_query Query filtered to include foreign samples & only selected swiss samples
downsample_swiss_sequences <- function (
  qcd_gisaid_query, db_connection, max_sampling_frac, favor_exposures = F, 
  verbose = T, subsample_by_canton = T, outdir = NULL
) {
  # get available sequences and bag exposure data
  qcd_gisaid_query_temp <- qcd_gisaid_query %>%
    filter(iso_country == "CHE") %>%
    mutate(week = date_trunc('week', date))
  bag_exposures <- get_bag_exposures(db_connection)
  
  # Get available seqs and cases numbers per week, possibly stratified by canton
  sampling_data_raw <- get_weekly_case_and_seq_data(
    db_connection = db_connection, 
    qcd_gisaid_query = qcd_gisaid_query, 
    by_canton = subsample_by_canton) 
  
  # restrict case data range to first and last week of available samples 
  # (this takes into account the analysis date range specified in qcd_gisaid_query)
  min_week <- qcd_gisaid_query_temp %>%
    summarize(min(week, na.rm = T)) %>%
    collect()
  max_week <- qcd_gisaid_query_temp %>%
    summarize(max(week, na.rm = T)) %>%
    collect()
    
  # Calculate # seqs per week
  # and divide proportionally amongst cantons if data stratified by canton
  sampling_data <- sampling_data_raw %>%
    mutate(week = as.character(week)) %>%
    filter(
      week <= max_week[[1]], week >= min_week[[1]], 
      (is.na(canton_code) | canton_code != "FL")) %>%  # don't count samples from FL, do count samples not associated with a canton
    group_by(week) %>%
    mutate(
      max_seqs_from_week = floor(sum(n_conf_cases) * max_sampling_frac),
      n_ideal_sample = quota_largest_remainder(
        votes = n_conf_cases, 
        n_seats = max_seqs_from_week[1]),
      n_to_sample = case_when(
        is.na(canton_code) ~ n_ideal_sample,  # hope that there are always enough Swiss-wide sequences to fill up the quota when the confirmed cases aren't attributed to a specific canton. If not, code will just sample all available sequences.
        T ~ pmin(n_ideal_sample, n_seqs_total))) %>%
    arrange(week, canton_code)
  
  # sample sequences 
  all_samples <- qcd_gisaid_query_temp %>%
    select(strain, week, division, gisaid_epi_isl, iso_country_exposure) %>%
    collect()
  sampled_strains <- c()
  for (i in 1:nrow(sampling_data)) {
    week_i <- sampling_data[[i, "week"]]
    division_i <- sampling_data[[i, "division"]]
    canton_code_i <- sampling_data[[i, "canton_code"]]
    n_samples_i <- as.numeric(sampling_data[[i, "n_to_sample"]])
    
    if (n_samples_i == 0) {  # if no samples desired, continue
      next
    }
    
    if (is.na(canton_code_i)) {  # if not stratifying by canton, or the canton is not known, take randomly from all Switzerland
      all_samples_i <- all_samples %>%
        filter(week == week_i, !(strain %in% sampled_strains))
    } else if (is.na(division_i)) {  # if no sequences available, continue
      if (verbose) {
        print(paste(
          "sampling", n_samples_i, 
          "samples out of 0", 
          "available from", case_when(is.na(canton_code_i) ~ "all Switzerland", T ~ canton_code_i),
          "in week", format(as.Date(week_i), "%Y-%m-%d")))
      }
      next
    } else {  # else take from the right canton
      all_samples_i <- all_samples %>%
        filter(week == week_i, division == division_i, !(strain %in% sampled_strains))
    }
    if (verbose) {
      print(paste(
        "sampling", n_samples_i, 
        "samples out of", nrow(all_samples_i), 
        "available from", ifelse(test = is.na(canton_code_i), yes = "all Switzerland", no = canton_code_i),
        "in week", format(as.Date(week_i), "%Y-%m-%d")))
    }
    n_samples_i <- min(n_samples_i, nrow(all_samples_i))  # correct for case when samples not attributed to a canton and there aren't enough sequences
    sampled_strains_i <- sample_strains(bag_exposures, all_samples_i, favor_exposures, n_samples_i)
    sampled_strains <- c(sampled_strains, sampled_strains_i)
  }

  # output downsampling data table and figure
  if (!(is.null(outdir))) {
    report_downsampling(sampling_data, outdir, max_sampling_frac)
  }
  
  # update master seq query to exclude non-sampled swiss strains
  qcd_gisaid_query <- qcd_gisaid_query %>%
    filter(iso_country != 'CHE' | strain %in% !! sampled_strains)
  return(qcd_gisaid_query)
}

sample_strains <- function(bag_exposures, all_samples_i, favor_exposures, n_samples_i) {
  if (favor_exposures) {
    # add bag meldeformular exposure data (not on GISAID) if there is any for samples from the week
    all_samples_i <- merge(
      x = bag_exposures, y = all_samples_i,
      by = "gisaid_epi_isl", 
      all.y = T) %>%
      mutate(iso_country_exposure = dplyr::coalesce(iso_country_exposure.x, iso_country_exposure.y)) 
    # shuffle samples but put foreign exposure samples first
    exp_samples_i <- all_samples_i %>%
      filter(!is.na(iso_country_exposure) & iso_country_exposure != "CHE")
    other_samples_i <- all_samples_i %>% 
      filter(is.na(iso_country_exposure) | iso_country_exposure == "CHE")
    shuffled_samples_i <- rbind(
      exp_samples_i[sample(nrow(exp_samples_i), replace = F), ],
      other_samples_i[sample(nrow(other_samples_i), replace = F), ]
    )
  } else {
    # shuffle samples
    shuffled_samples_i <- all_samples_i[sample(nrow(all_samples_i), replace = F), ]
  }
  sampled_strains_i <- shuffled_samples_i$strain[1:n_samples_i]
  return(sampled_strains_i)
}

report_downsampling <- function(sampling_data, outdir, max_sampling_frac) {
  write.csv(
    x = sampling_data, 
    file = paste(outdir, "swiss_downsampling_data.csv", sep = "/"),
    row.names = F)
  sampling_plot <- ggplot(
    data = sampling_data,
    aes(x = as.Date(week))) + 
    geom_col(aes(y = n_to_sample, fill = "Number of sequences analyzed")) + 
    geom_col(aes(y = -n_conf_cases * max_sampling_frac, 
                 fill = paste(max_sampling_frac * 100, "% of confirmed cases", sep = ""))) + 
    facet_wrap(canton_code ~ ., scales = "free_y") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b. %y") + 
    scale_y_continuous(limits = symmetric_limits) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_blank()) + 
    labs(x = element_blank(), y = "Count")
  ggsave(
    filename = paste(outdir, "swiss_downsampling.png", sep = "/"),
    plot = sampling_plot)
}

#' Get GISAID_EPI_ISL and exposure country for sequences with BAG-recoreded foreign
#' exposures. To be merged into GISAID data.
get_bag_exposures <- function(db_connection) {
  exposure_sql <- "select gisaid_epi_isl, iso_country_exp
      from gisaid_sequence gs
  join sequence_identifier si on si.gisaid_id = gs.gisaid_epi_isl
  left join viollier_test vt on si.ethid = vt.ethid
  left join bag_meldeformular bm on vt.sample_number = bm.sample_number
  where iso_country_exp is not null and
  iso_country_exp not in ('CHE', 'XXX');"
  res <- DBI::dbSendQuery(conn = db_connection, statement = exposure_sql)
  bag_exposures <- DBI::dbFetch(res = res)
  DBI::dbClearResult(res = res)
  return(bag_exposures %>% rename(iso_country_exposure = iso_country_exp))
}

# ------------------------------------------------------------------------------

#' Get number of cases by exposure country and month of case confirmation from 
#' BAG meldeformular.
#' Null entries for 'exp_land' or entries of 'Schweiz' in BAG meldeformular are not reported.
get_exposures_per_country_month <- function(db_connection, min_date, max_date) {
  print("Getting number of cases by exposure country and month of case confirmation from BAG meldeformular.")
  exposures_per_country_month <- dplyr::tbl(
    db_connection, "bag_meldeformular") %>%
    filter(!is.na(iso_country_exp), !(iso_country_exp %in% c('CHE', 'XXX')), 
           fall_dt <= !! max_date, 
           fall_dt >= !! min_date) %>%
    select(iso_country_exp, fall_dt) %>%
    collect() %>%
    mutate(date = format(fall_dt, "%Y-%m-01")) %>%
    group_by(date, iso_country_exp) %>%
    summarise(n_exposures = n()) %>%
    rename("iso_country" = "iso_country_exp")
  return(exposures_per_country_month)
}

#' Estimate number of tourists arriving in Switzerland by origin country and 
#' month of arrival. 
#' @param db_connection
#' @param max_date Character date, e.g. "2021-01-30" to extrapolate data to.
get_tourist_arrivals_per_country_month <- function(
  db_connection, min_date, max_date
) {
  print("Estimating number of tourists arriving in Switzerland by origin country and month of arrival.")
  tourist_arrivals_per_country_month <- dplyr::tbl(
    db_connection, "ext_fso_tourist_accommodation") %>%
    filter(date <= !! max_date, date >= !! min_date) %>%
    collect() %>%
    ungroup() %>%
    tidyr::complete(
      date = seq.Date(from = min(date), to = as.Date(max_date), by = "month"), 
      iso_country) %>%  
    group_by(iso_country) %>%
    arrange(date) %>%
    tidyr::fill(n_arrivals) %>% # fill in missing months with prev month's value
    rename(n_tourist_arrivals = n_arrivals)
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
    ungroup() %>%
    tidyr::complete(
      date = seq.Date(from = as.Date(min_date), to = as.Date(max_date), by = "month"), 
      iso_country, wirtschaftsabteilung) %>% 
    group_by(iso_country, wirtschaftsabteilung) %>%
    arrange(date) %>%
    tidyr::fill(n_permits) %>% # fill in months with quarterly values, missing quarters with previous quarter's value
    group_by(iso_country, date) %>%
    summarise(n_commuter_permits = sum(n_permits)) # sum across sectors
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
    select(iso_country, date, new_cases_per_million) %>%
    collect() %>%
    group_by(iso_country) %>%
    arrange(date) %>%
    mutate(cumul_cases_per_million = cumsum(new_cases_per_million)) %>%
    mutate(n_infectious_cases_per_million = lead(cumul_cases_per_million, n = 10) - cumul_cases_per_million) %>%
    mutate(date = as.Date(format(date, "%Y-%m-01"))) %>%
    filter(!is.na(n_infectious_cases_per_million)) %>%  # filters out NAs at ends of lead date ranges
    group_by(date, iso_country) %>%
    summarise(avg_daily_n_infectious_per_million = mean(n_infectious_cases_per_million)) %>% 
    ungroup() %>%
    tidyr::complete(
      date = seq.Date(from = min(date), to = as.Date(max_date), by = "month"), 
      iso_country) %>%  
    group_by(iso_country) %>%
    arrange(date) %>%
    tidyr::fill(avg_daily_n_infectious_per_million)  # fill in months with previous monthly value
  return(avg_infectious_per_country_month)  
}

#' Get data to estimate number of infectious arrivals into Switzerland by origin 
#' country and month of arrival.
#' @param db_connection
#' @param max_date Character date, e.g. "2021-01-30" to extrapolate data to.
get_infected_arrivals_per_country_month <- function(
  db_connection, min_date, max_date
) {
  print("Estimating number of infectious arrivals into Switzerland by origin country and month of arrival.")
  tourists <- get_tourist_arrivals_per_country_month(
    db_connection, min_date, max_date)
  commuters <- get_commuter_permits_per_country_month(
    db_connection, min_date, max_date)
  infectious <- get_avg_infectious_per_country_month(
    db_connection, min_date, max_date)
  arrivals <- merge(x = tourists, y = commuters, all = T, 
                    by = c("date", "iso_country")) %>%
    tidyr::replace_na(replace = list("n_tourist_arrivals" = 0, 
                                     "n_commuter_permits" = 0))
  infectious_arrivals <- merge(x = arrivals, y = infectious, 
                               by = c("date", "iso_country"), all.x = T)
  warn_missing_infectious_arrivals(infectious_arrivals)
  infectious_arrivals <- infectious_arrivals %>% 
    filter(!is.na(avg_daily_n_infectious_per_million)) 
  return(infectious_arrivals)
}

warn_missing_infectious_arrivals <- function(infectious_arrivals) {
  missing_case_data <- infectious_arrivals %>%
    filter(is.na(avg_daily_n_infectious_per_million), 
           date >= as.Date("2020-03-01")) %>%
    group_by(iso_country, date) %>%
    summarize(missing_case_data = T) 
  if (nrow(missing_case_data) > 0) {
    warning(
      "These countries are missing case data from after 2020-03-01 but have tourist/commuter arrivals into Switzerland. These arrivals will be ignored for context sequence subsampling.",
      paste(capture.output(print(missing_case_data)), collapse = "\n"))
  }
}

plot_estimated_travel_cases <- function(travel_cases, outdir) {
  print("Plotting estimated infectious arrivals per country and month.")
  require(ggplot2)
  top_10_countries <- travel_cases %>%
    group_by(iso_country) %>%
    summarize(total_travel_cases = sum(n_travel_cases)) %>%
    top_n(n = 8, wt = total_travel_cases) %>%
    select(iso_country)
  p <- ggplot(
    data = travel_cases %>% 
      mutate(source_country = case_when(
        iso_country %in% top_10_countries$iso_country ~ iso_country,
        T ~ "Other")),
    aes(x = format(date, "%Y-%m"), y = n_travel_cases)) +
    geom_col(aes(fill = source_country)) + 
    theme_bw() + 
    labs(x = element_blank(), y = "No. individuals") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  ggsave(
    filename = paste(outdir, "estimated_travel_cases.png", sep = "/"), 
    plot = p)
}

psum <- function(..., na.rm = T) { 
  rowSums(do.call(cbind,list(...)), na.rm = na.rm) 
} 

#' Get a set of strains proportional to estimated travel connections to Switzerland.
#' Estimated infectious individuals arriving and travel exposure cases are both
#' in absolute numbers, so sum for total number infectious inidiviuals arriving
#' in Switzerland from each country.
#' @param db_connection
#' @param travel_data_weights character with comma-delimited factors to scale
#' number of BAG-recorded exposures, tourist arrivals, and commuter permits by, 
#' respectively. Ex: '1,1,1'
#' @return Dataframe with information on number of infected individuals arriving
#' each month from each geographic origin. NA is only filled with zeros in col
#' n_travel_cases, otherwise NA means no source data.
get_travel_cases <- function(
  db_connection, min_date, max_date, outdir = NULL, travel_data_weights
) {
  print("Pulling tourist accommodation and cross-border commuter stats from FSO into the database.")
  source("database/R/import_fso_travel_stats.R")
  import_fso_tourist_accommodation(db_connection = db_connection)
  import_fso_cross_border_commuters(db_connection = db_connection)
  
  print("Querying tourist and exposure information from database.")
  infected_arrivals <- get_infected_arrivals_per_country_month(
    db_connection, min_date, max_date) 
  exposures <- get_exposures_per_country_month(db_connection, min_date, max_date)

  print(paste(
    "Parsing travel context prior weights for exposures, estimated infected arrivals:", 
    travel_data_weights))
  travel_data_weights <- as.numeric(strsplit(travel_data_weights, split = ",")[[1]])
  
  print("Merging arrival and exposure data, weighting estimates based on different data sources.")
  travel_cases <- merge(
    x = infected_arrivals, y = exposures, 
    all = T, by = c("date", "iso_country")) %>%
    mutate(
      n_unweighted_arrivals = psum(n_tourist_arrivals, n_commuter_permits),
      n_arrivals = psum(
        n_tourist_arrivals * travel_data_weights[2],
        n_commuter_permits * travel_data_weights[3]),
      n_infectious_arrivals = avg_daily_n_infectious_per_million * n_arrivals / 1E6,
      n_travel_cases = psum(
        n_exposures * travel_data_weights[1],
        n_infectious_arrivals,
        na.rm = T)) %>%
    arrange(date, desc(n_infectious_arrivals))
  
  to_remove <- travel_cases %>%
    filter(is.na(iso_country)) %>%
    select(iso_country, country, n_arrivals, n_exposures)
  if (length(to_remove) > 0) {
    warning("Not adding travel context sequences for these entries because no valid iso country code found.\n",
            paste(capture.output(print(to_remove)), collapse = "\n"))
    travel_cases <- travel_cases %>% filter(!is.na(iso_country))
  }
  
  print("Completing travel case data with zero values for months with no travel cases.")
  travel_cases <- travel_cases %>%
    select(-c(country, n_arrivals)) %>%
    ungroup() %>%
    tidyr::complete(
      date = seq.Date(from = min(date), to = as.Date(max_date), by = "month"), 
      iso_country) %>%  
    group_by(iso_country) %>%
    arrange(date) %>%
    tidyr::replace_na(replace = list(
      "n_travel_cases" = 0)) %>%  # fill in missing months with 0 value
    mutate(
      n_travel_cases = case_when(
        n_travel_cases < 0 ~ 0,
        T ~ n_travel_cases))  # when negative (can happen due to corrections in case count data), assume zero travel cases
  
  if (!is.null(outdir)) {
    plot_estimated_travel_cases(
      travel_cases = travel_cases,
      outdir = outdir)
    
    print(paste("Writing out results to", outdir))
    write.csv(
      x = travel_cases, 
      file = paste(outdir, "estimated_travel_cases.csv", sep = "/"), 
      row.names = F)
  }

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
    stop("No travel context sequences selected. Either increase timeframe or travel context scale factor.")
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
  if (verbose) {
    print(priorities_command)
  }
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
  python_path, reference, verbose = F
) {
  n_lineages <- nrow(lineages)
  initiated_df <- F
  for (i in 1:n_lineages) {
    lineage <- lineages$pangolin_lineage[i]
    lineages_included_str <- lineages$lineages_aggregated[i]
    lineages_included <- strsplit(lineages_included_str, split = ", ")[[1]]
    
    focal_seqs <- qcd_gisaid_query %>%
      filter(iso_country == "CHE", 
             pangolin_lineage %in% !! lineages_included) %>%
      select(strain) %>%
      collect()
    prospective_context_seqs <- qcd_gisaid_query %>%
      filter(iso_country != "CHE",
             pangolin_lineage %in% !! lineages_included) %>%
      select(strain, gisaid_epi_isl, iso_country, date_str, date) %>%
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
      reference = reference, 
      verbose = verbose)
    
    n_similarity_seqs <- ceiling(length(focal_strains) * similarity_context_scale_factor)
    similarity_strains_i <- priorities %>%
        arrange(desc(priority)) %>%
        mutate(priority_idx = 1:n()) %>%
        top_n(n = min(n_similarity_seqs, n()), wt = priority_idx) %>%
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
  initiated_df <- F
  for (i in 1:n_lineages) {
    lineage <- lineages$pangolin_lineage[i]
    lineages_included_str <- lineages$lineages_aggregated[i]
    lineages_included <- strsplit(lineages_included_str, split = ", ")[[1]]
    
    focal_strains_i <- qcd_gisaid_query %>%
      filter(iso_country == "CHE", 
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
      select(gisaid_epi_isl, strain, date_str, date, region, iso_country, division, 
             iso_country_exposure, nextstrain_clade, pangolin_lineage, 
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
    if (!initiated_df) {
      alignment_info <- alignment_info_i
      initiated_df <- T
    } else {
      alignment_info <- rbind(alignment_info, alignment_info_i)
    }
  }
  return(alignment_info %>% arrange(desc(n_seqs)))
}
  