# This script is to test whether the rate of new introductions from a country
# is effected by whether travellers from the country are under quarantine orders.

source("database/R/utility.R")
source("utility_functions.R")
source("generate_figures/functions.R")

s <- F
workdir <- "/Users/nadeaus/NonRepoProjects/cov-swiss-phylogenetics/grapevine/jan-dec_-01_max_sampling_-5_context-sf"
outdir <- paste(workdir, "output", sep = "/")
system(command = paste("mkdir -p", outdir))
db_connection <- open_database_connection()

#' TODO(?) extend time period to start of epidemic rather than start on first date 
#' of travel quarantine order and include indicator variable for CH border closure
#' (Mar 13 - Jun 15)
#' Date range for analysis is constrained to quarantine list data start/end dates.
test_travel_quarantine_effect <- function(
  s, workdir, outdir, db_connection, min_date, max_date
) {
  # indep variables: source country, whether travellers under quarantine order by day
  travel_quarantine_by_country_day <- get_travel_quarantine_by_country_day(
    db_connection = db_connection)  
  # indep variable: source country incidence per capita per day
  country_incidence <- get_country_incidence(
    min_date = min(model_data$date), 
    max_date = max(model_data$date),
    db_connection = db_connection)  

  warn_missing_incidence_data(travel_quarantine_by_country_day, country_incidence)
  travel_quarantine_by_country_day <- travel_quarantine_by_country_day %>%
    filter(iso_code %in% country_incidence$iso_code)  # countries on quarantine list without incidence data are dropped from analysis
  
  # dep variable: number of FOPH recorded exposures from each country each day
  recorded_exposures <- get_exposures_per_country_day(
    min_date = min(travel_quarantine_by_country_day$date), 
    max_date = max(travel_quarantine_by_country_day$date),
    db_connection = db_connection)
  # dep variable: number of lineages departing a country each day
  chains_asr <- load_chains_asr(s = F, workdir = workdir)
  chains_asr_representative <- pick_origin_representative_chains(chains_asr)
  chains_asr_long <- pivot_chains_longer(chains_asr_representative)
  
  indep_vars <- merge(
    x = travel_quarantine_by_country_day, y = country_incidence,
    by = c("iso_code", "date"), 
    all = T) 
  
  model_data_daily <- merge(
    x = indep_vars, y = recorded_exposures, 
    by = c("iso_code", "date"),
    all = T) %>%
    tidyr::replace_na(list(
      "quarantine_order" = F,  # any country not explicitly on the quarantine list has no quarantine order
      "n_exposures" = 0)) %>%  # if no exposure data recorded in bag_meldeformular for that country/date, assume 0 exposures
    group_by(iso_code) %>%
    mutate(quarantine_y_and_n = length(unique(quarantine_order)) == 2) %>%
    filter(quarantine_y_and_n) %>%  # keep only countries that were on the quarantine list at some point and off it at some point
    select(-quarantine_y_and_n)
  model_data_monthly <- model_data_daily %>%
    mutate(date_month = format(date, "%y-%m-01")) %>%
    group_by(date_month, iso_code) %>%
    summarize(
      n_exposures = sum(n_exposures),
      new_cases_per_million = sum(new_cases_per_million),
      quarantine_order = sum(quarantine_order) / n()) 
  
  # Model: daily # exposures ~ source country + source country incidence + quarantine indicator
  # Assumption: travel is constant year-round
  fitted_model <- fit_exposures_model(model_data = model_data_daily)  # A TRUE quarantine order on a day is associated with 0.03 more daily exposures across all countries, taking into account each country's baseline and independent of source country incidence
  fitted_model <- fit_exposures_model(model_data = model_data_monthly)  # Being on the quarantine list in a month is associated with 0.8 more monthly exposures across all countries, taking into account each country's baseline and independent of source country incidence
  fit_exposures_model(model_data_daily %>% filter(iso_code == "XKX"))
  
  # Report analysis results
  plot_travel_exposures_by_quarantine_status(
    travel_quarantine_by_country_day = travel_quarantine_by_country_day, 
    recorded_exposures = recorded_exposures, 
    outdir = outdir)
}

fit_exposures_model <- function(model_data, return_fitted_model = F) {
  # geographic distance as predictor
  # getting rid of country-specific variable gives more DOF
  fitted_model <- lm(
    data = model_data,
    formula = n_exposures ~ new_cases_per_million + quarantine_order)
  print(summary(fitted_model))
  if (return_fitted_model) {
    return(fitted_model)
  }
}



plot_travel_exposures_by_quarantine_status <- function(
  travel_quarantine_by_country_day, recorded_exposures, outdir
) {
  travel_exposures_by_quarantine_status <- ggplot(
    mapping = aes(x = date,
        y = iso_code)) + 
    geom_point(
      data = travel_quarantine_by_country_day,
      aes(color = quarantine_order)) + 
    geom_point(
      data = recorded_exposures,
      aes(size = n_exposures)) + 
    scale_color_manual(
      values = c(`FALSE` = "grey", `TRUE` = "red"),
      name = "Country under FOPH quarantine order") + 
    scale_size_continuous(
      name = "Number recorded exposures")
  
  ggsave(
    filename = paste(outdir, "travel_exposures_by_quarantine_status.png", sep = "/"),
    plot = travel_exposures_by_quarantine_status)
}



