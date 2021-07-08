#' Get a color assignment for countries so that colors are standardized across
#' plots. Assigns a color to all countries from database table gisaid_sequence.
#' Inspired by stack overflow: https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
#' @param db_connection
get_country_colors <- function(
  db_connection
) {
  require(RColorBrewer)
  countries <- dplyr::tbl(db_connection, "gisaid_sequence") %>%
    group_by(iso_country) %>%
    summarize(n_seqs = n()) %>%
    arrange(desc(n_seqs)) %>%
    collect()
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]  # I like these colors the best, they get priority
  all_colors = grDevices::colors()[grep(
    '(gr(a|e)y)|(black)|(white)|(seashell4)|(mediumorchid)', 
    grDevices::colors(), 
    invert = T)]  # These are a zillion colors, but some are really obnoxious 
  all_colors <- gplots::col2hex(all_colors)
  excluded_colors <- c("#BEAED4", "#CAB2D6", "#B3B3B3", "#FFFF99", "#E7298A")  # these are baddies - they look too similar to other colors
  
  set.seed(seed = 2600)
  
  initial_colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  initial_colors <- initial_colors[!(initial_colors %in% excluded_colors)]
  other_colors <- all_colors[!(all_colors %in% c(initial_colors, excluded_colors))]
  remaining_colors <- sample(
    other_colors, 
    size = (nrow(countries) - length(initial_colors)))
  
  colors <- c(initial_colors, remaining_colors)
  names(colors) <- countries$iso_country
  colors[['XXX']] <- "grey"
  colors[['OTHER']] <- "brown"
  colors[['UNKNOWN']] <- "white"

  return(colors)
}

#' Load grapevine results: one row per sample, columns describe transmission 
#' chain membership, transmission chain origin, transmission chain timeline.
#' @param workdir Directory with grapevine output.
#' @param min_chain_size Only return data for chains of >= min_chain_size.
#' @param viollier_only Only return data for samples from Viollier.
#' @return chains: dataframe with one row per chain; samples: dataframe with 
#' one row per focal sample; origins: dataframe with one row per node that
#' is a foreign attachment point of a chain. 
load_grapevine_results <- function(
  workdir, min_chain_size = 1, viollier_only = F
) {
  print(paste("Loading transmission chain data from", workdir))
  chains_with_asr <- rbind(
    load_chain_asr_data(l = F, workdir = workdir) %>%
      filter(size >= min_chain_size) %>%
      mutate(
        chain_idx = 1:n(),
        chains_assumption = "max",
        tree = gsub(tree, pattern = "_l_F|_l_T", replacement = "")),
    load_chain_asr_data(l = T, workdir = workdir) %>%
      filter(size >= min_chain_size) %>%
      mutate(
        chain_idx = 1:n(),
        chains_assumption = "min",
        tree = gsub(tree, pattern = "_l_F|_l_T", replacement = ""))
  ) 
  
  print("Formatting transmission chain data into per-sample, per-chain, and per-node dataframes.")
  # one row per chain
  chains_data <- chains_with_asr %>% select(
    tree, chains_assumption, chain_idx, 
    foreign_mrca, foreign_tmrca, foreign_tmrca_CI, 
    mrca, tmrca, tmrca_CI,
    size, tips, tip_nodes,
    n_intl_subclades
  )
  
  # one row per node that is a foreign_mrca of a focal chain
  origin_data <- chains_with_asr %>% 
    group_by(tree, chains_assumption, foreign_mrca, foreign_tmrca, foreign_tmrca_CI) %>% 
    mutate(n_chains_descending = n()) %>%
    group_by(n_chains_descending, .add = T) %>%
    summarize_at(
      .vars = vars(ends_with("_loc_weight")),
      .funs = list(~head(., 1))) %>%   # if multiple chains have same foreign mrca, e.g. at a polytomy, the asr data is duplicated
    ungroup()
  
  # one row per sample
  sample_metadata <- load_sample_metadata(workdir = workdir)
  samples_data <- pivot_chains_to_samples(
    chains = chains_data, metadata = sample_metadata) 
  if(!("iso_country" %in% colnames(samples_data))) {
    samples_data <- samples_data %>%
      mutate(
        iso_country = country_name_to_iso_code(country),
        iso_country_exposure = country_name_to_iso_code(country_exposure)
      )
  }
  samples_data <- samples_data %>% select(
    sample, tree, chains_assumption, chain_idx, sample_idx, gisaid_epi_isl, 
    strain, date_str, date, region, iso_country, division, iso_country_exposure, 
    nextstrain_clade, pangolin_lineage, originating_lab, submitting_lab, authors
  )
  
  if (viollier_only) {
    # take only samples from Viollier and submitted by us
    samples_data <- samples_data %>% filter(
      originating_lab == "Viollier AG",
      submitting_lab == "Department of Biosystems Science and Engineering, ETH Zürich")
    # prune chains to only samples from Viollier and submitted by us
    chains_data <- pivot_chains_to_samples(
      chains = chains_data, metadata = sample_metadata
    ) %>% filter(
      originating_lab == "Viollier AG",
      submitting_lab == "Department of Biosystems Science and Engineering, ETH Zürich"
    ) %>% group_by(
      tree, chains_assumption, chain_idx, 
      foreign_mrca, foreign_tmrca, foreign_tmrca_CI, 
      mrca, tmrca, tmrca_CI
    ) %>% summarize(
      size = n(),
      tips = paste0(sample, collapse = ", ")
    )
  }
  
  grapevine_results <- list(
    chains = chains_data,
    samples = samples_data,
    origins = origin_data
  )
  return(grapevine_results)
}


#' Join ancestral state location estimate for foreign mrca to focal chain data.
#' @param l Boolean indicating whether focal descendents of polytomies are the
#' same transmission chain.
#' @param chains_only If true, only returns chains data, without ASR.
#' @param chains_dirname Use to specify directory containing chains data per lineage.
#' @return Dataframe with fields describing transmission chains and normalized 
#' ancestral state location weights.
load_chain_asr_data <- function(
  l, workdir, chains_only = F, chains_dirname = "chains"
) {
  l_suffix <- paste("_l_", ifelse(test = l, yes = "T", no = "F"), sep = "")
  chains_suffix <- paste(l_suffix, "_chains.txt$", sep = "")
  print(paste("Looking for files ending in", chains_suffix, "within", paste(workdir, "tmp", chains_dirname, sep = "/")))
  chains_files <- list.files(
    path = paste(workdir, "tmp", chains_dirname, sep = "/"), pattern = chains_suffix)
  prefixes <- gsub(x = chains_files, pattern = chains_suffix, replacement = "")
  is_first <- T
  for (i in 1:length(prefixes)) {
    prefix <- prefixes[i]
    chains_file <- paste(workdir, "tmp", chains_dirname, chains_files[i], sep = "/")
    chains <- read.delim(file = chains_file, stringsAsFactors = F)
    
    if (!chains_only) {
      asr_filename <- paste(prefix, l_suffix, "_tree_data_with_asr.txt", sep = "")
      asr_file <- paste(workdir, "tmp/asr", asr_filename, sep = "/")
      asr <- read.delim(file = asr_file, stringsAsFactors = F, quote = "")
      chains_with_asr <- merge(
        x = chains, y = asr %>% select(node, ends_with("_loc_weight")),
        all.x = T, by.x = "foreign_mrca", by.y = "node") %>%
        mutate(tree = paste(prefix, l_suffix, sep = ""))
    } else {
      chains_with_asr <- chains %>%
        mutate(tree = paste(prefix, l_suffix, sep = ""))
    }
    
    if (is_first) {
      is_first <- F
      chains_all <- chains_with_asr
    } else {
      chains_all <- merge(x = chains_all, y = chains_with_asr, all = T)
    }
  }
  
  # Apply manual correction: LSD returns '2021' for date of tips with date '2020-12-31'
  chains_all[
    chains_all$size == 1 & 
      chains_all$tmrca == "2021" &
      grepl(x = chains_all$tips, pattern = "2020-12-31"), "tmrca"] <- "2020-12-31"
  chains_with_imprecise_mrca <- chains_all %>%
    filter(is.na(as.Date(chains_all$tmrca)))
  if (nrow(chains_with_imprecise_mrca) > 0) {
    warning(paste(
      "Some chains have imprecise tmrca estimate.",
      "LSD mysteriously replaces dates for tips on 2020-12-31 with 2021,", 
      "but this problem is already manually fixed in this function.", 
      "What else could be wrong?"))
    print(chains_with_imprecise_mrca)
  }
  return(chains_all)
}

#' @param workdir Directory containing grapevine's results
#' @return Dataframe with columns 'date', 'iso_code', 'ind_type' 
#' (tourist_arrivals or commuter_permits), and 'n_infectious_inds' (estimated 
#' number of infectious individuals of type ind_type arriving in Switzerland
#' in month of date).
load_origin_prior <- function(workdir) {
  origin_estimates <- read.csv(
    file = paste(
      workdir, "tmp/alignments/estimated_travel_cases.csv", 
      sep = "/"), 
    header = T, stringsAsFactors = F)
  
  origin_estimates_long <- origin_estimates %>%
    tidyr::pivot_longer(
      cols = c("n_tourist_arrivals", "n_commuter_permits", "n_exposures"),
      names_to = "ind_type",
      names_prefix = "n_",
      values_to = "n_inds") %>%
    mutate(n_infectious_inds = case_when(
      ind_type == "exposures" ~ 
        n_inds,
      ind_type %in% c("tourist_arrivals", "commuter_permits") ~ 
        avg_daily_n_infectious_per_million * n_inds / 1E6)) %>%
    select(date, iso_country, ind_type, n_infectious_inds) %>%
    filter(!is.na(n_infectious_inds), n_infectious_inds > 0)  # remove months and countries with no estimated infectious arrivals
  return(origin_estimates_long)
}

#' @param workdir Directory containing grapevine's results
#' @return Dataframe with columns 'date', 'iso_country', 'n_inds'.
load_travel_context <- function(workdir) {
  travel_context <- read.csv(
    file = paste(
      workdir, "tmp/alignments/travel_strains.txt", 
      sep = "/"), 
    header = T, stringsAsFactors = F)
  
  travel_context_summary <- travel_context %>%
    mutate(date = format(as.Date(date), "%Y-%m-01")) %>%
    group_by(date, iso_country) %>%
    summarize(n_inds = n())
  return(travel_context_summary)
}

#' Pivot chains dataframe to one row per sample with chain_idx information
#' Optionally attach sample metadata to table.
pivot_chains_to_samples <- function(
  chains, metadata = NULL, metadata_samplecol = "tree_label"
) {
  samples <- chains %>%
    tidyr::separate(
      col = tips,
      sep = ", ",
      into = paste("sample", 1:max(chains$size), sep = "_"),
      remove = T,
      fill = "right") %>% 
    tidyr::pivot_longer(
      cols = paste("sample", 1:max(chains$size), sep = "_"),
      values_to = "sample",
      names_to = "sample_idx",
      names_prefix = "sample_") %>%
    filter(!(is.na(sample)))
  if (!is.null(metadata)) {
    samples<- merge(
      x = samples, y = metadata,
      by.x = "sample", by.y = metadata_samplecol, all.x = T)
  }
  return(samples)
}

#' Load metadata from each alignment, concatenate
load_sample_metadata <- function(workdir, pattern = "*_metadata.tsv") {
  metadata_path <- paste(workdir, "tmp/alignments", sep = "/")
  metadata_files <- list.files(path = metadata_path, pattern = pattern, full.names = T)
  metadata <-
    metadata_files %>% 
    purrr::map_df(~readr::read_tsv(
      ., 
      quote = "\'",
      col_types = readr::cols(
        date = readr::col_date(format = "%Y-%m-%d")
      )
    ))
  print(paste("Loaded and concatenated", length(metadata_files), "metadata files."))
  
  # Remove duplicate entries for outgroup
  metadata_duplicates <- metadata %>% 
    group_by(gisaid_epi_isl) %>%
    summarize(n_occurances = n()) %>%
    filter(n_occurances > 1)
  for (i in 1:nrow(metadata_duplicates)) {
    epi_isl <- metadata_duplicates[i, "gisaid_epi_isl"]
    if (metadata_duplicates[i, "n_occurances"] != length(metadata_files)) {
      stop(paste("Strain", epi_isl, "found in multiple trees and is not in outgroup!"))
    } else {
      print(paste("Removing duplicate metadata entries for outgroup strain", epi_isl))
    }
  }
  metadata <- metadata %>% filter(!duplicated(metadata, fromLast = F))
  return(metadata)
}

#' Pivot ASR data from wide to long format, remove ASR assignments of Switzerland
#' for chains descending from a polytomy which have a the polytomy node as both 
#' mrca and foreign_mrca, replace periods in country names with spaces.
#' @param origins Dataframe with fileds describing normalized ancestral state 
#' location weights at nodes from which focal transmission chains descend.
#' @return Dataframe in long format with one column for chain origin and another 
#' for asr_contribution (normalized ancestral state location weights).
pivot_origins_longer <- function(origins) {
  origins_long <- origins %>% 
    mutate(origin_idx = 1:n()) %>%
    tidyr::pivot_longer(
    cols = ends_with("_loc_weight"),
    names_to = "origin",
    values_to = "asr_contribution") %>%
    mutate(
      origin = gsub(x = origin, pattern = "_loc_weight", replacement = "")) %>%
    # filter(origin != "Switzerland") %>%  # for l = T, a polytomy node can be both mrca and foreign_mrca. Count these chains as unclassifiable.
    mutate(origin = gsub(x = origin, pattern = "\\.", replacement = " "))
  return(origins_long)
}

#' @return Dataframe with n most common origins by total asr 
#' contribution across all chains for each chains assumption.
get_most_common_origins <- function(origins_long, n = 8) {
  most_common_origins <- origins_long %>%
    group_by(origin, chains_assumption) %>% 
    summarise(asr_contribution = sum(asr_contribution, na.rm = T)) %>%
    ungroup() %>% group_by(chains_assumption) %>%
    arrange(desc(asr_contribution)) %>%
    mutate(contribution_idx = 1:n()) %>%
    top_n(n, wt = -contribution_idx)
  return(most_common_origins %>% select(origin, chains_assumption))
}

#' Plot samples by inferred transmission chain
plot_chains <- function(
  workdir, outdir, country_colors, min_chain_size = 2, plot_height_in = 15
) {
  foreign_mrca_color <- "grey"
  mrca_color <- "black"
  
  chains <- load_chain_asr_data(l = F, workdir = workdir) %>%  # load max chains, then group by polytomy in plot
    filter(size > min_chain_size) %>%
    mutate(chain_idx = 1:n())
    
  sample_metadata <- load_sample_metadata(workdir = workdir)
  samples <- pivot_chains_to_samples(chains = chains, metadata = sample_metadata)
  
  samples_info <- samples %>%
    group_by(chain_idx) %>%
    mutate(first_sample_in_chain = min(date),
           last_sample_in_chain = max(date)) %>%
    ungroup() %>%
    tidyr::unite(col = "unique_foreign_mrca", tree_pangolin_lineage, foreign_mrca)  # add facetting variable for each separate foreign attachment node on global tree
  
  chain_info <- samples_info %>% 
    group_by(chain_idx, first_sample_in_chain, last_sample_in_chain, unique_foreign_mrca, tmrca) %>%
    summarize(size = n()) %>%
    full_join(y = chains) %>%
    group_by(tree, foreign_mrca) %>%
    arrange(desc(tmrca)) %>%
    mutate(chain_order_within_foreign_mrca = 1:n()) %>%
    ungroup()
  
  samples_info <- samples_info %>%
    left_join(y = chain_info) %>%  # attach chain_order_within_foreign_mrca info to samples
    mutate(date = as.Date(date))  # format date information for samples
  
  chain_info <- chain_info %>%
    tidyr::separate(
      foreign_tmrca_CI, 
      into = c("foreign_tmrca_CI_min", "foreign_tmrca_CI_max"),
      sep = ", ") %>%
    mutate(
      foreign_tmrca_CI_min = gsub(x = foreign_tmrca_CI_min, pattern = "c\\(", replacement = ""),
      foreign_tmrca_CI_max = gsub(x = foreign_tmrca_CI_max, pattern = "\\)", replacement = ""),
      foreign_tmrca_CI_min = as.Date(foreign_tmrca_CI_min),
      foreign_tmrca_CI_max = as.Date(foreign_tmrca_CI_max)) %>%
    tidyr::separate(
      tmrca_CI,
      into = c("tmrca_CI_min", "tmrca_CI_max"),
      sep = ", ") %>%
    mutate(
      tmrca_CI_min = gsub(x = tmrca_CI_min, pattern = "c\\(", replacement = ""),
      tmrca_CI_max = gsub(x = tmrca_CI_max, pattern = "\\)", replacement = ""),
      tmrca_CI_min = as.Date(tmrca_CI_min),
      tmrca_CI_max = as.Date(tmrca_CI_max),
      first_sample_in_chain = as.Date(first_sample_in_chain),
      last_sample_in_chain = as.Date(last_sample_in_chain),
      foreign_tmrca = as.Date(foreign_tmrca),
      tmrca = as.Date(tmrca))  # format date information for chains
  
  samples_info <- samples_info %>% 
    mutate(unique_foreign_mrca = factor(
      x = samples_info$unique_foreign_mrca,
      levels = unique(unlist(samples_info %>% arrange(foreign_tmrca) %>% select(unique_foreign_mrca)))))  # order facets by node times
  chain_info <- chain_info %>% 
    mutate(unique_foreign_mrca = factor(
      x = chain_info$unique_foreign_mrca,
      levels = unique(unlist(chain_info %>% arrange(foreign_tmrca) %>% select(unique_foreign_mrca)))))  # order facets by node times
  
  custom_theme_elements <- theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(), 
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    panel.spacing.x = unit(0,"line"),
    strip.text.y = element_blank(),
    legend.position = c(0.07, 0.45),
    legend.background = element_rect(fill = "transparent"),
    panel.spacing = unit(0, "lines"))
  axis_labs <- labs(x = element_blank(), y = element_blank())
  
  lockdown_end <- as.Date("2020-06-15")
  lockdown_start <- as.Date("2020-03-25")

  transmission_chain_plot <- ggplot() + 
    geom_segment(
      data = chain_info,
      aes(x = first_sample_in_chain, 
          xend = last_sample_in_chain,
          y = chain_order_within_foreign_mrca, 
          yend = chain_order_within_foreign_mrca),
      linetype = "dashed") +
    geom_point(
      data = chain_info,
      aes(x = foreign_tmrca,
          y = chain_order_within_foreign_mrca),
      shape = 4, color = foreign_mrca_color) +
    geom_point(
      data = chain_info,
      aes(x = tmrca,
          y = chain_order_within_foreign_mrca),
      shape = 4, color = mrca_color) +
    geom_errorbarh(
      data = chain_info,
      aes(xmin = tmrca_CI_min,
          xmax = tmrca_CI_max,
          y = chain_order_within_foreign_mrca),
      color = mrca_color) +
    geom_errorbarh(
      data = chain_info,
      aes(xmin = foreign_tmrca_CI_min,
          xmax = foreign_tmrca_CI_max,
          y = chain_order_within_foreign_mrca),
      color = foreign_mrca_color) +
    geom_point(
      data = samples_info,
      aes(x = date,
          y = chain_order_within_foreign_mrca)) +
          # color = country_exposure)) +
    # # Plot colored points over the top so they stand out more
    # geom_point(
    #   data = samples_info %>% filter(country_exposure != country),
    #   aes(x = date, 
    #       y = chain_order_within_foreign_mrca, 
    #       color = country_exposure)) + 
    # geom_rect(
    #   data = chain_info,
    #   aes(xmin = lockdown_start, 
    #       xmax = lockdown_end,
    #       ymin = min(chain_order_within_foreign_mrca),
    #       ymax = max(chain_order_within_foreign_mrca)),
    #   alpha = 0.2) + 
    geom_vline(xintercept = lockdown_start, linetype = "longdash") + 
    geom_vline(xintercept = lockdown_end, linetype = "longdash") + 
    facet_grid(unique_foreign_mrca ~ ., scales = "free_y", space = "free_y") + 
    custom_theme_elements + 
    scale_x_date(date_breaks = "month", date_labels = "%b %y") + 
    axis_labs

  pdf(NULL)
  ggsave(transmission_chain_plot, 
         file = paste(outdir, "transmission_chains.png", sep = "/"), 
         height = plot_height_in)
  dev.off()
  
  # + geom_scatterpie(
  #     data = mrca_asr_data_2,
  #     aes(x = x_axis_days, y = cluster_order_midpoint, r = 4),
  #     cols = c(countries_to_plot, "other"),
  #     color = NA) +
  #   theme_bw() + 
  #   custom_theme_elements +
  #   scale_fill_manual(
  #     values = color_scale, 
  #     name = element_blank(), 
  #     aesthetics = c("fill", "color"),
  #     breaks = c(countries_to_plot[countries_to_plot != "Switzerland"], "other")) + 
  #   x_scale +
  #   labs(y = "\n\n", x = element_blank()) + 
  #   facet_grid(foreign_mrca ~ ., scales = "free_y", space = "free_y")
}

#' Plot estimated origins of focal transmission lineages through time (prior vs. 'posterior').
#' @param grapevine_results Pre-loaded results of grapevine pipeline.
#' @param country_colors Named list with color values for countries.
plot_chain_origins <- function(
  workdir, grapevine_results, outdir = NULL, country_colors, origins_to_color = NULL
) {
  # Count each origin once per unique foreign MRCA
  origins_data <- get_origin_prior_vs_posterior_data(
    workdir = workdir,
    origins_representative = grapevine_results$origins)
  
  # Count each origin once per unique chain
  origins_per_chain <- merge(
    x = grapevine_results$origins,
    y = grapevine_results$chains,
    by = c("tree",  "chains_assumption", "foreign_mrca", "foreign_tmrca", "foreign_tmrca_CI"),
    all = T)
  origins_data_per_chain <- get_origin_prior_vs_posterior_data(
    workdir = workdir,
    origins_representative = origins_per_chain)
  
  # Count each origin once per unique chain, use the focal mrca data and not the foreign mrca date
  origins_data_per_chain_mrca <- get_origin_prior_vs_posterior_data(
    workdir = workdir,
    origins_representative = origins_per_chain,
    posterior_date = "tmrca")
  
  if (is.null(origins_to_color)) {
    most_common_origins <- get_most_common_origins_to_plot(
      origins_data = origins_data)
    most_common_origins_per_chain <- get_most_common_origins_to_plot(
      origins_data = origins_data_per_chain)
    origins_to_color = most_common_origins$iso_country
    origins_per_chain_to_color = most_common_origins_per_chain$iso_country
  } else {
    origins_to_color <- origins_to_color
    origins_per_chain_to_color <- origins_to_color
  }
  
  origins_plot <- make_origins_plot(
    country_colors = country_colors, 
    origins_data = origins_data, 
    origins_to_color = origins_to_color)
  origins_per_chain_plot <- make_origins_plot(
    country_colors = country_colors, 
    origins_data = origins_data_per_chain, 
    origins_to_color = origins_per_chain_to_color)
  origins_per_chain_tmrca_plot <- make_origins_plot(
    country_colors = country_colors, 
    origins_data = origins_data_per_chain_mrca,
    origins_to_color = origins_per_chain_to_color)
  
  if (!is.null(outdir)) {
    save_origins_plot(
      plot = origins_plot,
      outdir = outdir,
      plotname = "chain_origins_per_foreign_mrca_foreign_tmrca.png")
    save_origins_plot(
      plot = origins_per_chain_plot,
      outdir = outdir,
      plotname = "chain_origins_per_chain_foreign_tmrca.png")
    save_origins_plot(
      plot = origins_per_chain_tmrca_plot,
      outdir = outdir,
      plotname = "chain_origins_per_chain_tmrca.png")
  } else {
    return(list("chain_origins_per_foreign_mrca_foreign_tmrca" = origins_plot,
                "chain_origins_per_chain_foreign_tmrca" = origins_per_chain_plot,
                "chain_origins_per_chain_tmrca" = origins_per_chain_tmrca_plot))
  }
}

#' @return Long-format dataframe with prior, and posterior under different 
#' chains assumptions, estimates for transmission chain sources.
get_origin_prior_vs_posterior_data <- function(
  workdir, origins_representative, posterior_date = "foreign_tmrca"
) {
  # Load data
  origin_prior_long <- load_origin_prior(workdir = workdir) 
  travel_context <- load_travel_context(workdir = workdir)
  origins_posterior_long <- pivot_origins_longer(origins_representative)
  
  # Fill in Unknown for origins that are unclassifiable (under max clusters assumption where mrca == foreign_mrca or when a chain clusters only with similarity context sequences)
  origins_summary <- origins_posterior_long %>% 
    group_by(chains_assumption, tree, foreign_mrca, !!sym(posterior_date), origin_idx) %>% 
    summarize(total_asr_contribution = sum(asr_contribution, na.rm = T))
  
  unknown_origins <- origins_summary[origins_summary$total_asr_contribution == 0, ] %>%
    mutate(origin = "UNKNOWN", asr_contribution = 1) %>%
    select(-c(total_asr_contribution))
  origins_posterior_long_2 <- merge(
    x = origins_posterior_long, y = unknown_origins, all = T, 
    by = c("chains_assumption", "tree", "foreign_mrca", posterior_date, "origin_idx", "asr_contribution", "origin")) %>%
    mutate(origin = recode(origin, "CHE" = "UNKNOWN"))
  
  origins_summary_2 <- origins_posterior_long_2 %>% 
    group_by(chains_assumption, tree, foreign_mrca, origin_idx) %>% 
    summarize(total_asr_contribution = sum(asr_contribution, na.rm = T))  # now all origins have asr_contributions totaling 1
  origins_posterior_long <- origins_posterior_long_2
  
  # Define factors so that complete will complete data for all levels of country, etc. (necessary for area plot)
  all_months <- sort(unique(c(
    format(as.Date(origins_posterior_long[[posterior_date]]), "%Y-%m-01"),
    origin_prior_long$date)))
  prior_countries <- unique(origin_prior_long$iso_country)
  origin_prior_long$date <- factor(
    x = origin_prior_long$date,
    levels = all_months)
  travel_context$date <- factor(
    x = travel_context$date,
    levels = all_months)
  origins_posterior_long_countries_not_in_prior <- origins_posterior_long %>%
    filter(!(origin %in% c(prior_countries, "UNKNOWN")), asr_contribution > 0)
  if (nrow(origins_posterior_long_countries_not_in_prior) > 0) {
    print(origins_posterior_long_countries_not_in_prior)
    stop("Some countries have asr_contribution even though they are not in the prior!")
  }
  origins_posterior_long <- origins_posterior_long %>% 
    filter(origin %in% c(prior_countries, "UNKNOWN")) %>%
    mutate(
      date = format(as.Date(!!sym(posterior_date)), "%Y-%m-01"),
      date = factor(x = date, levels = all_months),
      iso_country = factor(x = origin, levels = c(prior_countries, "UNKNOWN")))
  travel_context$iso_country <- factor(
    x = travel_context$iso_country,
    levels = c(prior_countries, "UNKNOWN"))
  
  # Merge data into long format
  origin_prior_vs_posterior_data <- merge(
    x = origin_prior_long %>%
      ungroup() %>%
      tidyr::complete(
        date, ind_type, iso_country,
        fill = list(n_infectious_inds = 0)
      ) %>%
      mutate(
        date_type = "estimate_month",
        est_type = "prior_base_estimates") %>%
      rename(
        "n_inds" = "n_infectious_inds",
        "n_inds_type" = "ind_type"),
    y = origins_posterior_long %>%
      group_by(date, iso_country, chains_assumption) %>%
      summarize(asr_contribution = sum(asr_contribution, na.rm = T)) %>%
      ungroup() %>%
      tidyr::complete(
        date, iso_country, chains_assumption,
        fill = list(asr_contribution = 0)
      ) %>%
      mutate(
        date_type = paste(posterior_date, "month", sep = "_"),
        n_inds_type = "sum_fractional_lineage_asr",
        est_type = "posterior") %>%
      rename(
        "n_inds" = "asr_contribution") %>%
      select(
        chains_assumption, iso_country, date, date_type, n_inds, n_inds_type,
        est_type),
    all = T,
    by = c("iso_country", "date", "date_type", "n_inds", "n_inds_type", "est_type")
  )
  origin_prior_vs_posterior_data_w_travel_context <- merge(
    x = origin_prior_vs_posterior_data,
    y = travel_context %>%
      ungroup() %>%
      tidyr::complete(
        date, iso_country,
        fill = list(n_inds = 0)
      ) %>% mutate(
        date_type = "sampling_month",
        n_inds_type = "n_seqs",
        est_type = "prior_sequences_used"),
    all = T,
    by = c("iso_country", "date", "date_type", "n_inds", "n_inds_type", "est_type")
  )
  return(origin_prior_vs_posterior_data_w_travel_context %>% ungroup())
}

get_most_common_origins_to_plot <- function(
  origin_freq_threshold = 0.10, origin_count_threshold = 3, origins_data
) {
  origins_freq <- origins_data %>%
    group_by(est_type, chains_assumption, date) %>%
    mutate(total_n_inds = sum(n_inds)) %>%
    ungroup() %>%
    group_by(est_type, chains_assumption, date, iso_country) %>%
    summarize(freq_inds = sum(n_inds) / total_n_inds[[1]]) %>% # NaN means 0/0, so n_total_inds 0
    ungroup()
  
  origins_count <- origins_data %>%
    group_by(est_type, chains_assumption, iso_country) %>%
    summarize(total_n_inds = sum(n_inds))
  
  most_freq_origins <- origins_freq %>%
    filter(freq_inds > origin_freq_threshold)
  most_count_origins <- origins_count %>%
    filter(total_n_inds > origin_count_threshold)
  most_common_origins <- most_freq_origins %>%
    filter(iso_country %in% most_count_origins$iso_country)
  
  cat(
    length(unique(most_common_origins$iso_country)), 
    "origins will be highlighted in the plot for acheiving >",
    origin_freq_threshold * 100,
    "% frequency in some month in the prior or either posterior and at least",
    origin_count_threshold, "total lineages:\n",
    paste0(unique(most_common_origins$iso_country), collapse = ", "), 
    "\n"
  )
  return(most_common_origins)
}

make_origins_plot <- function(origins_data, country_colors, origins_to_color) {
  origins_data_to_plot <- origins_data %>%
    mutate(
      origin_to_plot = case_when(
        iso_country %in% origins_to_color ~ iso_country,
        T ~ "OTHER"
      )
    ) %>%
    tidyr::unite(
      col = "group", origin_to_plot, n_inds_type, remove = F
    ) %>% tidyr::unite(
      col = "facet", est_type, chains_assumption, remove = F
    ) %>% group_by(
      date, origin_to_plot, group, facet, n_inds_type) %>%
    summarize(n_inds = sum(n_inds)) 
  
  origins_data_to_plot$n_inds_type <- factor(
    x = origins_data_to_plot$n_inds_type,
    levels = c("sum_fractional_lineage_asr", "n_seqs", "exposures", "tourist_arrivals", "commuter_permits")
  )
  origins_data_to_plot$facet <- factor(
    x = origins_data_to_plot$facet,
    levels = c(
      "prior_base_estimates_NA", 
      "prior_sequences_used_NA",
      "posterior_max", 
      "posterior_min"),
    labels = c(
      "Prior base estimates", 
      "Travel context\n sequences used",
      "Tree-based estimates\n (Smallest plausible chains)",
      "Tree-based estimates\n (Largest plausible chains)"))
  
  plot <- ggplot(
    data = origins_data_to_plot,
    aes(x = as.Date(date), y = n_inds)) + 
    facet_grid(
      facet ~ .,
      scales = "free") + 
    geom_area(
      aes(group = group, fill = as.factor(origin_to_plot), alpha = n_inds_type),
      position = "stack") +  # 'fill' means height corresponds to percentage of all classifiable foreign mrcas in month
    scale_fill_manual(
      values = country_colors,
      labels = iso_code_to_country_name,
      name = "Origin country") +
    scale_alpha_manual(
      values = c(
        sum_fractional_lineage_asr = 1,
        n_seqs = 1,
        exposures = 0.75, 
        tourist_arrivals = 0.5, 
        commuter_permits = 0.25),
      labels = c(
        sum_fractional_lineage_asr = "Number of lineages (estimated)",
        n_seqs = "Number of sequences",
        exposures = "Reported exposure",
        tourist_arrivals = "Infected tourist (estimated)",
        commuter_permits = "Infected commuter (estimated)"),
      name = "Import type") + 
    labs(
      x = element_blank(),
      y = "Number of imports") + 
    scale_x_date(date_labels = "%b. %Y", date_breaks = "1 month") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  return(plot)
}

save_origins_plot <- function(
  plotname, plot, outdir = "figures", single_col_width = 11.4) {
  pdf(NULL)
  ggsave(
    plot = plot + 
      theme(
        legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin()) + 
      guides(
        fill = guide_legend(
          ncol = 3, byrow = T, title.position = "top", title.hjust = 0.5),
        alpha = guide_legend(
          nrow = 2, byrows = T, title.position = "top", title.hjust = 0.5)),
    file = paste(outdir, plotname, sep = "/"),
    width = single_col_width,
    height = single_col_width * 2,
    units = "cm"
  )
  dev.off()
}

#' Table prior vs. posterior transmission chain origins in different time periods.
#' @param origins One row per foreign mrca of focal transmission chains,
#' with a column giving how many chains descend. 
table_chain_origins <- function(
  workdir, outdir = NULL, period_breaks = c("2020-05-01"), 
  n_origins_to_table = 15, origins, return_table = F
) {
  # Load prior 
  prior_origins_long <- load_travel_context(workdir = workdir)
  # Load posterior
  posterior_origins_long <- pivot_origins_longer(origins)
  
  # Summarize prior and posterior origin contributions by period
  dates_range <- range(c(prior_origins_long$date, posterior_origins_long$foreign_tmrca))
  complete_period_breaks <- unique(c(
    min(dates_range, period_breaks),
    max(dates_range, period_breaks),
    period_breaks))  # most extreme dates + specified intermediate breakpoints
  posterior_by_period <- posterior_origins_long %>%
    mutate(
      period = cut.Date(x = as.Date(foreign_tmrca), 
                   breaks = as.Date(complete_period_breaks),
                   include.lowest = T)) %>%
    group_by(period, origin) %>%
    summarise(asr_contribution = sum(asr_contribution, na.rm = T)) %>%
    group_by(period) %>%
    mutate(frac_posterior = asr_contribution / sum(asr_contribution)) %>%
    rename(iso_country = origin)
  prior_by_period <- prior_origins_long %>%
    mutate(
      period = cut.Date(x = as.Date(date), 
                        breaks = as.Date(complete_period_breaks),
                        include.lowest = T)) %>%
    group_by(period, iso_country) %>%
    summarise(n_inds = sum(n_inds, na.rm = T)) %>%
    group_by(period) %>%
    mutate(frac_prior = n_inds / sum(n_inds))
  posterior_vs_prior_by_period <- merge(
    x = posterior_by_period, y = prior_by_period, by = c("period", "iso_country"),
    all = T) %>%
    mutate(frac_posterior_to_frac_prior_ratio = frac_posterior / frac_prior)
  
  # Format results table for top posterior origins, write out 
  most_common_origins <- get_most_common_origins(
    origins_long = posterior_origins_long, n = n_origins_to_table)
  posterior_vs_prior_by_period_to_table <- posterior_vs_prior_by_period %>%
    filter(iso_country %in% most_common_origins$origin) %>%
    select(period, frac_posterior_to_frac_prior_ratio, iso_country) %>%
    mutate(
      frac_posterior_to_frac_prior_ratio = signif(
        frac_posterior_to_frac_prior_ratio, digits = 2),
      country = iso_code_to_country_name(iso_code = iso_country)) %>%
    tidyr::pivot_wider(
      names_from = period,
      values_from = frac_posterior_to_frac_prior_ratio,
      names_prefix = "Period beginning ")
  
  if (!is.null(outdir)) {
    filename <- paste("posterior_vs_prior_origins_l_",
                      ifelse(test = l, yes = "T", no = "F"), ".csv", sep = "")
    write.csv(
      x = posterior_vs_prior_by_period_to_table,
      file = paste(outdir, filename, sep = "/"),
      row.names = F)
  }
  if (return_table) {
    return(posterior_vs_prior_by_period_to_table)
  }
}

#' Plot a lineage annotated with BAG meldeformular data. This function is an unfinished TODO
plot_lineage_with_exposure_location <- function(db_connection, lineage, workdir) {
  sql <- "select bm.exp_land, si.gisaid_id from
  sequence_identifier si
  left join viollier_test vt on si.ethid = vt.ethid
  left join bag_meldeformular bm on vt.sample_number = bm.sample_number
  where si.gisaid_id is not null"
  bag_exposures <- DBI::dbGetQuery(conn = db_connection, statement = sql)
  
  lineage <- "B.1.1.162"
  tree_filename <- paste(lineage, ".timetree.nex", sep = "")
  tree <- treeio::read.beast(file = paste(workdir, "tmp/lsd", tree_filename, sep = "/"))
  tree_data_with_asr_filename <- paste(lineage, "_m_3_p_1_l_F_tree_data_with_asr.txt", sep = "")
  tree_data_with_asr <- read.delim(file = paste(workdir, "tmp/asr", tree_data_with_asr_filename, sep = "/")) 
  
  asr_pie_data <- tree_data_with_asr %>% select(node, all_of(ends_with("loc_weight")))
  rownames(asr_pie_data) <- tree_data_with_asr$strain
  pies <- ggtree::nodepie(asr_pie_data, cols = 2:ncol(asr_pie_data))
  
  tree_2 <- full_join(tree, tree_data_with_asr)
  
  p <-  ggtree(tr = tree_2) +
    geom_inset(insets = pies) +
    geom_tiplab(aes(
      label = strain,
      color = country),
      hjust = -0.01)

    pdf(NULL)
    ggsave(
    file = paste(outdir, paste(prefix, "tree_with_asr.png", sep = "_"), sep = "/"), 
    plot = p)
    dev.off()
}

#' Plot barchart of sampling intensity through time.
#' TODO: include unsampled seqeunces from viollier again
#' @param dates_to_highlight Vector of character dates to draw vertical lines for.
plot_sampling_intensity <- function(
  db_connection, workdir, outdir, max_date, max_sampling_frac, 
  dates_to_highlight = NULL
) {
  sample_metadata <- load_sample_metadata(workdir = workdir)
  samples_query <- dplyr::tbl(db_connection, "gisaid_sequence") %>%
    filter(gisaid_epi_isl %in% !! sample_metadata$gisaid_epi_isl)
  weekly_case_and_seq_data <- get_weekly_case_and_seq_data(
    db_connection = db_connection, qcd_gisaid_query = samples_query) %>%
    filter(as.Date(week) <= as.Date(max_date))
  plot_data <- weekly_case_and_seq_data %>%
    mutate(n_unseq_conf_cases = pmax(n_conf_cases - n_seqs_viollier - n_seqs_other, 0)) %>%  # For any weeks where # sequences > # confirmed cases, set # unsequenced confirmed cases to 0.
    tidyr::pivot_longer(
      cols = c("n_seqs_viollier", "n_seqs_other", "n_unseq_conf_cases"),
      names_to = "case_type",
      values_to = "n_cases", 
      names_prefix = "is_|n_") %>%
    mutate(week = as.Date(week))
  plot_data$case_type <- factor(
    x = plot_data$case_type, 
    levels = c("unseq_conf_cases", "seqs_other", "seqs_viollier"))
  
  shared_theme <- theme_bw()
  
  freq_plot <- ggplot(
    data = plot_data, 
    aes(x = week, y = n_cases, fill = case_type)) + 
    geom_col(position = position_fill()) +
    labs(x = element_blank(), y = "Weekly percentage of cases")
  abs_plot <- ggplot(
    data = plot_data, 
    aes(x = week, y = n_cases, fill = case_type)) + 
    geom_col(position = position_stack()) +
    labs(x = element_blank(), y = "Weekly number of cases")

  add_vlines <- function(plot, dates_to_highlight) {
    if (is.null(dates_to_highlight)) {
      return(plot)
    } else {
      dates_to_highlight <- as.Date(dates_to_highlight, "%Y-%m-%d")
      for (date in dates_to_highlight) {
        plot <- plot + geom_vline(xintercept = date, linetype = "dashed")
      }
      return(plot)
    }
  }
  
  weekly_samples_vs_cases <- ggplot(
    data = weekly_case_and_seq_data, 
    aes(x = as.Date(week), y = n_seqs_total / n_conf_cases)) + 
    geom_col() + 
    geom_hline(yintercept = max_sampling_frac, color = "red", linetype = "dashed") + 
    geom_text(
      aes(
        y = case_when(
          (n_seqs_total / n_conf_cases) > 0.01 ~  (n_seqs_total / n_conf_cases) - 0.0001,
          T ~ (n_seqs_total / n_conf_cases) + 0.0001),
        hjust = case_when(
          (n_seqs_total / n_conf_cases) > 0.01 ~ 1,
          T ~ 0),
        color = case_when(
          (n_seqs_total / n_conf_cases) > 0.01 ~ "inside_color",
          T ~ "outside_color"),
        label = paste(n_seqs_total, "/", n_conf_cases)),
      angle = 90, vjust = 0.5) + 
    labs(x = "Week", y = "Empirical sampling fraction") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") + 
    scale_color_manual(values = c("inside_color" = "white", "outside_color" = "black")) + 
    shared_theme + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  weekly_samples_vs_cases <- weekly_samples_vs_cases %>% 
    add_vlines(dates_to_highlight = dates_to_highlight)

  shared_scale_fill <- scale_fill_manual(
      values = RColorBrewer::brewer.pal(n = 3, name = "Dark2"),
      labels = c("seqs_viollier" = "Sequenced, from Viollier", 
                 "seqs_other" = "Sequenced, from another lab",
                 "unseq_conf_cases" = "Unsequenced confirmed case"),
      name = element_blank())
  
  sampling_intensity_plot <- ggpubr::ggarrange(
    abs_plot + shared_scale_fill + shared_theme,
    freq_plot + shared_scale_fill + shared_theme,
    nrow = 2, common.legend = T, legend = "right")

  pdf(NULL)
  ggsave(
    file = paste(outdir, "sampling_intensity.png", sep = "/"), 
    plot = sampling_intensity_plot)
  ggsave(
    file = paste(outdir, "weekly_samples_vs_cases.png", sep = "/"),
    plot = weekly_samples_vs_cases)
  dev.off()
    
  return(plot_data)
}

#' Plot maximum and minimum-bound estimates on the number of newly sampled chains
#' per week (~introductions) and final chain samples per week (~extinctions).
#' @param last_sample_to_extinction_delay Days between the last sample in a 
#' transmission chain and the date at which the chain is recorded as having gone
#' extinct.
plot_introductions_and_extinctions <- function(
  workdir, outdir, min_chain_size = 1, last_sample_to_extinction_delay = 14
) {
  samples <- load_grapevine_results(
    workdir = workdir, min_chain_size = min_chain_size, viollier_only = T)$samples
  week_to_week_start <- data.frame(
    day = seq.Date(
    from = min(samples$date), 
    to = max(samples$date), 
    by = "day")) %>%
    mutate(week = format(day, "%Y.%W")) %>%
    group_by(week) %>%
    top_n(n = 1, wt = desc(day))
  
  first_samples <- samples %>%
    group_by(chains_assumption, chain_idx) %>%
    arrange(desc(date)) %>%
    mutate(sample_idx = 1:n()) %>%
    top_n(n = 1, wt = sample_idx) %>%
    mutate(event_type = "introduction")
  last_samples <- samples %>% 
    group_by(chains_assumption, chain_idx) %>%
    arrange(date) %>%
    mutate(sample_idx = 1:n()) %>%
    top_n(n = 1, wt = sample_idx) %>% 
    mutate(event_type = "extinction",
           date = date + last_sample_to_extinction_delay)
  introductions_and_extinctions <- rbind(first_samples, last_samples) %>%
    select(chains_assumption, date, chain_idx, event_type) %>%
    mutate(week = format(date, "%Y.%W")) %>%
    group_by(week, chains_assumption, event_type) %>%
    summarise(n_events = n()) %>%
    tidyr::pivot_wider(
      names_from = chains_assumption, 
      values_from = n_events,
      names_prefix = "events_",
      values_fill = list(n_events = 0)) %>%
    left_join(y = week_to_week_start, by = "week") %>%
    mutate(event_type = factor(event_type, levels = c("introduction", "extinction")))
    
  introductions_and_extinctions_plot <- ggplot() + 
    geom_errorbar(
      data = introductions_and_extinctions,
      aes(x = day, 
          ymin = events_min,
          ymax = events_max,
          color = event_type)) + 
    scale_x_date(date_breaks = "month", date_labels = "%b %y") + 
    scale_color_manual(
      values = c(extinction = "red", introduction = "blue"),
      labels = c(introduction = "New chain sampled",
                 extinction = "Chain presumed extinct"),
      name = "Event type") + 
    theme_bw() + 
    labs(x = element_blank(), y = "Weekly number of events")

  pdf(NULL)
  ggsave(
    file = paste(outdir, "introductions_and_extinctions.png", sep = "/"), 
    plot = introductions_and_extinctions_plot)
  dev.off()
  return(introductions_and_extinctions_plot)
}

#' @return Dataframe with columns date (at daily frequency), iso_code, 
#' quarantine_order (boolean), and n_regions_quarantine (NA for non-neighboring 
#' countries).
get_travel_quarantine_by_country_day <- function(db_connection) {
  foph_quarantine_list <- dplyr::tbl(
    db_connection, "foph_travel_quarantine") %>%
    collect()
  
  travel_quarantine_by_country_day <- foph_quarantine_list %>%
    rename("date" = "date_effective") %>%
    group_by(iso_code, date) %>%
    summarize(n_regions = n()) %>%
    ungroup() %>%
    mutate(
      n_regions_quarantine = case_when(
        iso_code %in% c("AUT", "FRA", "ITA", "DEU") ~ n_regions),
      quarantine_order = T,
      date = as.Date(date)) %>%
    select(iso_code, date, quarantine_order, n_regions_quarantine) %>%
    tidyr::complete(
      date, iso_code, 
      fill = list("quarantine_order" = F, "n_regions_quarantine" = NA)) %>%
    tidyr::complete(date = seq.Date(min(date), max(date), by = "day"), iso_code) %>%
    group_by(iso_code) %>%
    arrange(date) %>%
    tidyr::fill(quarantine_order)
  
  return(travel_quarantine_by_country_day)
}

#' @param iso_country List of country codes; if NULL (default) returns all.
get_country_incidence <- function(db_connection, min_date, max_date, iso_countries = NULL) {
  country_incidence_unfiltered <- dplyr::tbl(db_connection, "ext_owid_global_cases") %>%
    select(iso_country, date, new_cases_per_million, new_cases) %>%
    filter(date >= min_date, date <= max_date) %>%
    collect()
  if (is.null(iso_countries)) {
    country_incidence_filtered <- country_incidence_unfiltered
  } else {
    country_incidence_filtered <- country_incidence_unfiltered %>%
      filter(iso_country %in% !! iso_countries)
  }
  country_incidence <- country_incidence_filtered %>%
    tidyr::complete(date = seq.Date(min_date, max_date, by = "day"), iso_country) %>%
    mutate(
      new_cases_per_million = case_when(
        is.na(new_cases_per_million) ~ 0,  # assume days without data have 0 new cases
        new_cases_per_million < 0 ~ 0,  # replace days with negative incidence with 0
        T ~ new_cases_per_million),
      new_cases = case_when(
        is.na(new_cases) ~ 0,  # assume days without data have 0 new cases
        new_cases < 0 ~ 0,  # replace days with negative prevalence with 0
        T ~ as.numeric(new_cases)))  
  return(country_incidence)
}

warn_missing_incidence_data <- function(
  travel_quarantine_by_country_day, country_incidence
) {
  countries_missing_case_data <- unique(travel_quarantine_by_country_day$iso_code[!(
    travel_quarantine_by_country_day$iso_code %in% country_incidence$iso_code)])
  if (length(countries_missing_case_data) > 0) {
    missing_countries <- data.frame(
      iso_code = countries_missing_case_data,
      country_name = iso_code_to_country_name(countries_missing_case_data))
    warning("These countries are on the quarantine list but are missing from ",
            "database table ext_owid_global_cases and will therefore be excluded ",
            "from test travel quarantine effect analysis:\n",
            paste(capture.output(print(missing_countries)), collapse = "\n"))
  }
}

#' Get number of cases by exposure country and date of case confirmation from 
#' BAG meldeformular.
#' Null/unknown entries for 'iso_exp_land' or entries of 'CHE' in BAG meldeformular are not reported.
get_exposures_per_country_day <- function(db_connection, min_date, max_date) {
  print("Getting number of cases by exposure country and date of case confirmation from BAG meldeformular.")
  exposures_per_country_day <- dplyr::tbl(
    db_connection, "bag_meldeformular") %>%
    filter(!is.na(iso_country_exp), iso_country_exp != "CHE", iso_country_exp != "XXX", 
           fall_dt <= !! max_date, fall_dt >= !! min_date) %>%
    select(iso_country_exp, fall_dt) %>%
    collect() %>%
    rename("date" = "fall_dt") %>%
    group_by(date, iso_country_exp) %>%
    summarise(n_exposures = n())
  
  return(exposures_per_country_day)
}

#' Get distance from each country to Switzerland
#' @param non_focal_countries list of country codes to get distance from
#' @param focal_countries countries to get distance to
#' @param db_connection
#' @return dataframe of distances from each non-focal country centroid to 
#' the focal country centroids in km
get_country_distances <- function(
  non_focal_countries, focal_countries, db_connection
) {
  country_coordinates <- dplyr::tbl(db_connection, "ext_country_coordinates") %>%
    filter(iso_code %in% !! c(non_focal_countries, focal_countries)) %>%
    collect()
  all_countries <- c(non_focal_countries, focal_countries)
  missing_countries <- all_countries[!(all_countries %in% country_coordinates$iso_code)]
  if (length(missing_countries) > 1) {
    warning(paste("These countries are missing from geographic location dataset.", 
                  "They will not be included in distance output.", 
                  paste0(iso_code_to_country_name(missing_countries), collapse = ", ")))
  }
  non_focal_coordinates <- country_coordinates %>% 
    filter(iso_code %in% non_focal_countries)
  focal_coordinates <- country_coordinates %>% 
    filter(iso_code %in% focal_countries)
  distances_km <- as.data.frame(geodist::geodist(
    x = non_focal_coordinates,
    y = focal_coordinates,
    measure = "geodesic") / 1000)
  colnames(distances_km) <- paste("dist_to_", focal_coordinates$iso_code, sep = "")
  distances_km$iso_code <- non_focal_coordinates$iso_code
  return(distances_km)
}

#' Plot n_exposures or n_lineages by quarantine status. Used when testing for
#' quarantine effects on imports.
plot_dep_var_by_quarantine_status <- function(
  model_data, outdir = NULL, dep_varname, return_plot = T, plot_height_cm = 11.4
) {
  if (is.numeric(model_data$quarantine_order)) {
    model_data$quarantine_order <- model_data$quarantine_order >= 0.5
  }
  
  dep_var_by_quarantine_status <- ggplot(
    data = model_data,
    mapping = aes(
      x = date,
      y = iso_code)) + 
    geom_point(
      aes(color = quarantine_order)) + 
    geom_point(
      aes(
        size = ifelse(!!sym(dep_varname) == 0, NA, !!sym(dep_varname))),
      shape = 1) +
    scale_color_manual(
      values = c(`FALSE` = "grey", `TRUE` = "red"),
      name = "Country under FOPH quarantine order") +
    scale_size_continuous(
      name = dep_varname,
      range = c(0, 6)) + 
    theme_bw() + 
    labs(x = element_blank(), y = "Origin country")
  
  if (dep_varname == "n_lineages") {
    dep_var_by_quarantine_status <- dep_var_by_quarantine_status + 
      facet_grid(. ~ chains_assumption)
  }
  
  if (!is.null(outdir)) {
    outfile <- paste(dep_varname, "_by_quarantine_status.png", sep = "")
    pdf(NULL)
    ggsave(
      filename = paste(outdir, outfile, sep = "/"),
      plot = dep_var_by_quarantine_status,
      height = plot_height_cm, units = "cm")
    dev.off()
  }
  if (return_plot) {
    return(dep_var_by_quarantine_status)
  }
}

#' Plot the time-scaled phylogeny with tip locations, estimated ancestral node
#' locations, estimated focal transmission chains, and genetic similarity &
#' travel context tips.
#' @param tree Grapevine output, e.g. /Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/jan-dec_-01_max_sampling_-5_context-sf/tmp/lsd/B.1.509.timetree.nex
#' @param tree_data_with_asr Grapevine output, e.g. /Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/jan-dec_-01_max_sampling_-5_context-sf/tmp/asr/B.1.509_m_3_p_1_l_F_tree_data_with_asr.txt
#' @param chains Grapevine output, e.g. /Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/jan-dec_-01_max_sampling_-5_context-sf/tmp/chains/B.1.509_m_3_p_1_l_F_chains.txt
plot_tree <- function(
  tree, tree_data_with_asr, chains, country_colors, outdir = NULL, prefix = "tree"
) {
  
  country_colors[["CHE"]] <- "red"
  country_colors_for_pies <- country_colors
  names(country_colors_for_pies) <- paste(names(country_colors), "_loc_weight", sep = "")
  
  tree_data_to_plot <- tree_data_with_asr %>% 
    mutate(
      node_annotation = "",
      tip_type = "")
      # node_annotation = case_when(
      #   node %in% chains$mrca ~ "focal MRCA",
      #   T ~ ""),
      # tip_type  = case_when(
      #   travel_context ~ "travel context",
      #   similarity_context ~ "genetic context",
      #   focal_sequence ~ "",
      #   T ~ "root"))
  
  # Apply manual correction: LSD returns '2021' for date of tips with date '2020-12-31'
  tree_data_to_plot$date <- as.character(tree_data_to_plot$date)
  tree_data_to_plot[tree_data_to_plot$date == "2021" & !is.na(tree_data_to_plot$label) & grepl(x = tree_data_to_plot$label, pattern = "2020-12-31"), "date"] <- "2020-12-31"
  tree_data_to_plot <- tree_data_to_plot %>% filter(node > 0)
  
  tree_to_plot <- full_join(tree, tree_data_to_plot)
  
  loc_weight_cols <- which(grepl(x = colnames(tree_data_to_plot), pattern = "_loc_weight"))
  node_location_pies <- ggtree::nodepie(
    tree_data_to_plot, 
    cols = loc_weight_cols)
  node_location_pies_colored <- lapply(
    X = node_location_pies, 
    FUN = function(g) g + scale_fill_manual(values = country_colors_for_pies))
  
  p <- ggtree(
    tr = tree_to_plot,
    mrsd = max(as.Date(tree_data_to_plot$date)), 
    as.Date = T) +
    # geom_inset(insets = node_location_pies_colored) +
    geom_nodelab(aes(
      label = node_annotation)) +
    geom_tiplab(aes(
      label = iso_code_to_country_name(iso_country),
      # label = case_when(
      #   tip_type != "" ~ paste(iso_code_to_country_name(iso_country), tip_type, sep = ": "),
      #   T ~ iso_code_to_country_name(iso_country)),
      color = iso_country),
      hjust = -0.1, size = 2) +
    scale_color_manual(values = country_colors) +
    scale_x_date(
      date_labels = "%b. %Y",
      limits = c(as.Date("2019-12-01"), as.Date("2021-06-01"))) + 
    theme_tree2() + 
    theme(legend.position = "none")
    
  if (is.null(outdir)) {
    return(p)
  } else {
    pdf(NULL)
    ggsave(
      file = paste(outdir, paste(prefix, "_with_asr.pdf", sep = ""), sep = "/"),
      plot = p,
      dpi = 400
    )
    dev.off()
  }
}
