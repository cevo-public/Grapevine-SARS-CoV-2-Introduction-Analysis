#' Get a color assignment for countries so that colors are standardized across
#' plots. Assigns a color to top 30 countries from gisaid_sequence.
#' @param db_connection
#' @param name_type Country name convention to use for list names. One of 'english_name' or 'iso3'.
#' @param n_unique_colors [Optional] Number of countries to get unique colors. Others get grey.
get_country_colors <- function(
  db_connection, name_type = "english_name", n_unique_colors = NULL
) {
  countries <- dplyr::tbl(db_connection, "gisaid_sequence") %>%
    group_by(country) %>%
    summarize(n_seqs = n()) %>%
    arrange(desc(n_seqs)) %>%
    collect()
  countries <- countries %>% 
    mutate(iso3 = countrycode::countrycode(
      sourcevar = country, origin = "country.name", destination = "iso3c"))
  if (is.null(n_unique_colors)) {
    n <- nrow(countries)
    all_colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    set.seed(seed = 10)
    colors <- sample(all_colors, n)
    colors <- c(colors, "grey")
  } else {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_unique_colors)
    colors <- c(colors, rep("grey", 1 + nrow(countries) - n_unique_colors))
  }
  if (name_type == "english_name") {
    names(colors) <- c(countries$country, "other")
  } else if (name_type == "iso3") {
    names(colors) <- c(countries$iso3, "other")
  }
  return(colors)
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
load_sample_metadata <- function(workdir, pattern = "*_metadata.csv") {
  metadata_path <- paste(workdir, "tmp/alignments", sep = "/")
  metadata_files <- list.files(path = metadata_path, pattern = pattern, full.names = T)
  metadata <-
    metadata_files %>% 
    purrr::map_df(~readr::read_csv(
      ., 
      col_types = readr::cols(
        date = readr::col_date(format = "%Y-%m-%d"),
        date_str = readr::col_character()
      )
    ))
  print(paste("Loaded and concatenated", length(metadata_files), "metadata files."))
  return(metadata)
}

#' Plot samples by inferred transmission chain
plot_chains <- function(
  workdir, outdir, country_colors, min_chain_size = 2, plot_height_in = 15
) {
  foreign_mrca_color <- "grey"
  ch_mrca_color <- "black"
  
  chains <- load_chains_asr(s = F, workdir = workdir) %>%  # load max chains, then group by polytomy in plot
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
    group_by(chain_idx, first_sample_in_chain, last_sample_in_chain, unique_foreign_mrca, ch_tmrca) %>%
    summarize(size = n()) %>%
    full_join(y = chains) %>%
    group_by(tree, foreign_mrca) %>%
    arrange(desc(ch_tmrca)) %>%
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
      ch_tmrca_CI, 
      into = c("ch_tmrca_CI_min", "ch_tmrca_CI_max"),
      sep = ", ") %>%
    mutate(
      ch_tmrca_CI_min = gsub(x = ch_tmrca_CI_min, pattern = "c\\(", replacement = ""),
      ch_tmrca_CI_max = gsub(x = ch_tmrca_CI_max, pattern = "\\)", replacement = ""),
      ch_tmrca_CI_min = as.Date(ch_tmrca_CI_min),
      ch_tmrca_CI_max = as.Date(ch_tmrca_CI_max),
      first_sample_in_chain = as.Date(first_sample_in_chain),
      last_sample_in_chain = as.Date(last_sample_in_chain),
      foreign_tmrca = as.Date(foreign_tmrca),
      ch_tmrca = as.Date(ch_tmrca))  # format date information for chains
  
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
      aes(x = ch_tmrca,
          y = chain_order_within_foreign_mrca),
      shape = 4, color = ch_mrca_color) +
    geom_errorbarh(
      data = chain_info,
      aes(xmin = ch_tmrca_CI_min,
          xmax = ch_tmrca_CI_max,
          y = chain_order_within_foreign_mrca),
      color = ch_mrca_color) +
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
  
  ggsave(transmission_chain_plot, 
         file = paste(outdir, "transmission_chains.png", sep = "/"), 
         height = plot_height_in)
  
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

plot_source_prior <- function(
  workdir, outdir, country_colors, n_sources_to_plot = 10
) {
  source_estimates <- read.table(
    file = paste(
      workdir, "tmp/alignments/estimated_monthly_infectious_arrivals.txt", 
      sep = "/"), 
    sep = "\t", header = T, stringsAsFactors = F)
  
  source_estimates_long <- source_estimates %>%
    tidyr::pivot_longer(
      cols = c("n_tourist_arrivals", "n_commuter_permits"),
      names_to = "ind_type",
      names_prefix = "n_",
      values_to = "n_inds") %>%
   select(-c(n_infectious_arrivals, n_arrivals)) %>%
   mutate(n_infectious_inds = avg_daily_n_infectious_per_million * n_inds / 1E6)
  
  source_summary <- source_estimates_long %>%
    group_by(origin) %>% 
    summarise(prior_contribution = sum(n_infectious_inds, na.rm = T)) %>%
    arrange(desc(prior_contribution))  
  
  sources_to_plot <- source_estimates_long %>%
    mutate(
      source_to_plot = case_when(
        origin %in% source_summary$origin[1:n_sources_to_plot] ~ origin,
        T ~ "other")) %>%
    tidyr::unite(col = "source_to_plot_ind_type", source_to_plot, ind_type, remove = F)
  
  plot <- ggplot(
    data = sources_to_plot, 
    aes(x = as.Date(date), y = as.numeric(n_infectious_inds))) + 
    geom_area(
      aes(group = source_to_plot_ind_type, fill = source_to_plot, alpha = ind_type),
      position = "stack") +  # 'fill' means height corresponds to percentage of all classifiable foreign mrcas in month
    scale_fill_manual(values = country_colors, name = "Source country") + 
    scale_alpha_manual(
      values = c(tourist_arrivals = 1, commuter_permits = 0.65), 
      name = "Arrival type") + 
    labs(x = "Month of entry", 
         y = "Estimated number of infectious arrivals") + 
    scale_x_date(date_breaks = "month", date_labels = "%b %y") + 
    theme_bw()
  
  ggsave(plot = plot, file = paste(outdir, "chain_source_prior.png", sep = "/"))
}

#' Plot estimated origins of swiss transmission lineages through time.
plot_chain_sources <- function(s, workdir, outdir, country_colors) {
  chains_asr <- load_chains_asr(s = s, workdir = workdir)
  plot <- make_plot_chain_sources(chains_asr, country_colors)
  filename <- paste("chain_sources_s_", 
                    ifelse(test = s, yes = "T", no = "F"), ".png", 
                    sep = "")
  ggsave(plot = plot, file = paste(outdir, filename, sep = "/"))
}

#' Plot estimated origins of swiss transmission lineages through time.
#' @param tree_data_with_asr Tree data in tidytree treedata format with extra 
#' columns for ancestral state estimates.
#' @param chains Data frame with fields describing transmission chains.
#' @param country_colors Named list with color values for countries.
#' @return plot
#' TODO: fix that "other" shows up as transparent in plot
make_plot_chain_sources <- function(
  chains_asr, country_colors, n_sources_to_plot = 10
) {
  chains_asr_representative <- chains_asr %>% 
    mutate(foreign_mrca_month = format(as.Date(foreign_tmrca), "%Y-%m-01")) %>%
    group_by(tree, foreign_mrca) %>%
    mutate(child_idx = 1:n()) %>%
    top_n(n = 1, wt = child_idx)  # keep only one representative chain for chains descending from same polytomy so that polytomy ASR's aren't overweighted when s = F.
  
  chains_asr_long <- chains_asr_representative %>% tidyr::pivot_longer(
    cols = ends_with("_loc_weight"),
    names_to = "source",
    values_to = "asr_contribution") %>%
    mutate(
      source = gsub(x = source, pattern = "_loc_weight", replacement = "")) %>%
    filter(source != "Switzerland")  # for s = T, a polytomy node can be both ch_mrca and foreign_mrca. Count these chains as unclassifiable.
  
  source_summary <- chains_asr_long %>%
    group_by(source) %>% 
    summarise(asr_contribution = sum(asr_contribution, na.rm = T)) %>%
    arrange(desc(asr_contribution))  
  
  chains_asr_to_plot <- chains_asr_long %>%
    group_by(source, foreign_mrca_month) %>%
    summarise(asr_contribution = sum(asr_contribution, na.rm = T)) %>%
    mutate(
      source_to_plot = case_when(
        source %in% source_summary$source[1:n_sources_to_plot] ~ source,
        T ~ "other"),
      source_to_plot = gsub(x = source_to_plot, pattern = "\\.", replacement = " "))
  
  n_nodes_data <- chains_asr_representative %>% 
    group_by(foreign_mrca_month) %>%
    summarize(
      n_foreign_mrcas = n(), 
      n_foreign_mrcas_with_asr = sum(!is.na(Switzerland_loc_weight) & Switzerland_loc_weight != 1)) %>%  # Sum across nodes with NA for all locations, including Switzerland (no travel context sequences = no ASR) and nodes with Switzerland estimated as ASR (a polytomy where the same node is both ch_mrca and foreign_mrca).
    mutate(label = paste(n_foreign_mrcas_with_asr, "(", n_foreign_mrcas, ")", sep = ""))
  
  plot <- ggplot(
    data = chains_asr_to_plot, 
    aes(x = as.Date(foreign_mrca_month), y = asr_contribution)) + 
    geom_area(
      aes(group = source_to_plot, fill = source_to_plot),
      position = "stack") +  # 'fill' means height corresponds to percentage of all classifiable foreign mrcas in month
    geom_text(
      data = n_nodes_data,
      aes(y = 0, label = label),
      vjust = 1) +  # annotate bars with total number of foreign mrcas in month, including those unclassifiable for reasons commented above
    scale_fill_manual(values = country_colors) + 
    labs(x = "Month of foreign attachment point", 
         y = "Number of lineages") + 
    scale_x_date(date_breaks = "month", date_labels = "%b %y") + 
    theme_bw()

  return(plot)
}

#' Join ancestral state location estimate for foreign mrca to swiss chain data.
#' @param s Boolean indicating whether swiss descendents of polytomies are the 
#' same transmission chain.
#' @return Dataframe with fields describing transmission chains and normalized 
#' ancestral state location weights.
load_chains_asr <- function(s, workdir) {
  s_suffix <- paste("_s_", ifelse(test = s, yes = "T", no = "F"), sep = "")
  chains_suffix <- paste(s_suffix, "_chains.txt$", sep = "")
  chains_files <- list.files(
    path = paste(workdir, "tmp/chains", sep = "/"), pattern = chains_suffix)
  prefixes <- gsub(x = chains_files, pattern = chains_suffix, replacement = "")
  is_first <- T
  for (i in 1:length(prefixes)) {
    prefix <- prefixes[i]
    chains_file <- paste(workdir, "tmp/chains", chains_files[i], sep = "/")
    asr_filename <- paste(prefix, s_suffix, "_tree_data_with_asr.txt", sep = "")
    asr_file <- paste(workdir, "tmp/asr", asr_filename, sep = "/")
    asr <- read.delim(file = asr_file, stringsAsFactors = F)
    chains <- read.delim(file = chains_file, stringsAsFactors = F)
    chains_asr_tree <- merge(
      x = chains, y = asr %>% select(node, ends_with("_loc_weight")),
      all.x = T, by.x = "foreign_mrca", by.y = "node") %>%
      mutate(tree = paste(prefix, s_suffix, sep = ""))
    if (is_first) {
      is_first <- F
      chains_asr <- chains_asr_tree
    } else {
      chains_asr <- merge(x = chains_asr, y = chains_asr_tree, all = T)
    }
  }
  return(chains_asr)
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
  tree_data_with_asr_filename <- paste(lineage, "_m_3_p_1_s_F_tree_data_with_asr.txt", sep = "")
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

    ggsave(
    file = paste(outdir, paste(prefix, "tree_with_asr.png", sep = "_"), sep = "/"), 
    plot = p)
}

#' Plot barchart of sampling intensity through time.
#' TODO: include unsampled seqeunces from viollier again
plot_sampling_intensity <- function(db_connection, workdir, outdir, max_date) {
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

  weekly_samples_vs_cases <- ggplot(
    data = weekly_case_and_seq_data,
    aes(x = n_conf_cases, y = n_seqs_total)) + 
    geom_point(aes(color = week)) +
    labs(x = "Weekly number confirmed cases", y = "Weekly number sequences analyzed") + 
    shared_theme

  shared_scale_fill <- scale_fill_manual(
      values = RColorBrewer::brewer.pal(n = 3, name = "Dark2"),
      labels = c("seqs_viollier" = "Sequenced, from Viollier", 
                 "seqs_other" = "Sequenced, from another lab",
                 "unseq_conf_cases" = "Unsequenced confirmed case"),
      name = element_blank())
  
  sampling_intensity_plot <- ggpubr::ggarrange(
    freq_plot + shared_scale_fill + shared_theme,
    abs_plot + shared_scale_fill + shared_theme,
    nrow = 2, common.legend = T, legend = "right")
  
  ggsave(
    file = paste(outdir, "sampling_intensity.png", sep = "/"), 
    plot = sampling_intensity_plot)
  ggsave(
    file = paste(outdir, "weekly_samples_vs_cases.png", sep = "/"),
    plot = weekly_samples_vs_cases)
    
  return(plot_data)
}

#' Plot maximum and minimum-bound estimates on the number of introductions per
#' week.
plot_introductions <- function(workdir, outdir) {
  chains_max <- load_chains_asr(s = F, workdir = workdir) %>%  
    filter(size > min_chain_size) %>%
    mutate(chain_idx = 1:n())
  chains_min <- load_chains_asr(s = T, workdir = workdir) %>%  
    filter(size > min_chain_size) %>%
    mutate(chain_idx = 1:n())
  
  sample_metadata <- load_sample_metadata(workdir = workdir)
  samples <- rbind(
    pivot_chains_to_samples(chains = chains_max, metadata = sample_metadata) %>%
      mutate(chains_assumption = "max"),
    pivot_chains_to_samples(chains = chains_min, metadata = sample_metadata) %>%
      mutate(chains_assumption = "min")) %>%
    filter(originating_lab == "Viollier AG")
  
  week_to_week_start <- data.frame(
    day = seq.Date(
    from = min(samples$date), 
    to = max(samples$date), 
    by = "day")) %>%
    mutate(week = format(day, "%Y.%W")) %>%
    group_by(week) %>%
    top_n(n = 1, wt = desc(day))
  
  introduction_detections <- samples %>%
    group_by(chains_assumption, chain_idx) %>%
    arrange(desc(date)) %>%
    mutate(sample_idx = 1:n()) %>%
    top_n(n = 1, wt = sample_idx) %>%
    select(chains_assumption, date, chain_idx) %>%
    mutate(week = format(date, "%Y.%W")) %>%
    group_by(week, chains_assumption) %>%
    summarise(n_introductions = n()) %>%
    tidyr::pivot_wider(
      names_from = chains_assumption, 
      values_from = n_introductions,
      names_prefix = "introductions_",
      values_fill = list(n_introductions = 0)) %>%
    left_join(y = week_to_week_start, by = "week")
    
  introduction_plot <- ggplot() + 
    geom_errorbar(
      data = introduction_detections,
      aes(x = day, 
          ymin = introductions_min,
          ymax = introductions_max)) + 
    scale_x_date(date_breaks = "month", date_labels = "%b %y") + 
    theme_bw() + 
    labs(x = element_blank(), y = "Weekly number newly sampled introductions")
  
  ggsave(
    file = paste(outdir, "introductions.png", sep = "/"), 
    plot = introduction_plot)
}

