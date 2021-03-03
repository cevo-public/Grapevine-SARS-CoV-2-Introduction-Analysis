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
make_plot_chain_sources <- function(chains_asr, country_colors) {
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
        source %in% source_summary$source[1:8] ~ source,
        T ~ "other"))
  
  n_nodes_data <- chains_asr_representative %>% 
    group_by(foreign_mrca_month) %>%
    summarize(
      n_foreign_mrcas = n(), 
      n_foreign_mrcas_with_asr = sum(!is.na(Switzerland_loc_weight))) %>%  # nodes without ASR have NA for all locations, including Switzerland
    mutate(label = paste(n_foreign_mrcas_with_asr, "(", n_foreign_mrcas, ")", sep = ""))
  
  plot <- ggplot(
    data = chains_asr_to_plot, 
    aes(x = foreign_mrca_month, y = asr_contribution)) + 
    geom_col(aes(fill = source_to_plot), position = "fill") +  # bars show percentage of all foreign mrcas
    geom_text(
      data = n_nodes_data,
      aes(y = 1, label = label), 
      vjust = 0) +  # annotate bars with number of foreign mrcas
    scale_fill_manual(values = country_colors)
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
#' TODO: plot for samples actually included, rather than pulling eligible samples from the database
#' TODO: include unsampled seqeunces from viollier again
plot_sampling_intensity <- function(db_connection, qcd_gisaid_query, outdir) {
  weekly_case_and_seq_data <- get_weekly_case_and_seq_data(
    db_connection = db_connection, qcd_gisaid_query = qcd_gisaid_query) 
  plot_data <- weekly_case_and_seq_data %>%
    mutate(n_unseq_conf_cases = n_conf_cases - n_seqs_viollier - n_seqs_other) %>%
    tidyr::pivot_longer(
      cols = c("n_seqs_viollier", "n_seqs_other", "n_unseq_conf_cases"),
      names_to = "case_type",
      values_to = "n_cases", 
      names_prefix = "is_|n_") %>%
    mutate(week = as.Date(week), n_cases = as.numeric(n_cases), 
           n_conf_cases = as.numeric(n_conf_cases))
  
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
  abs_scatterplot <- ggplot(
    data = plot_data %>% filter(case_type == "seqs_viollier"),
    aes(x = n_conf_cases, y = n_cases)) + 
    geom_point(aes(color = week)) + 
    geom_abline(slope = 0.01)  # show what a specific percent sampling would look like
  show(abs_scatterplot)
  
  shared_scale_fill <- scale_fill_manual(
      values = RColorBrewer::brewer.pal(n = 3, name = "Dark2"),
      labels = c("seqs_viollier" = "Sequenced, from Viollier", 
                 "seqs_other" = "Sequenced, from another lab",
                 "unseq_conf_cases" = "Unsequenced confirmed case"),
      name = element_blank())
  shared_theme <- theme_bw()
  
  sampling_intensity_plot <- ggpubr::ggarrange(
    freq_plot + shared_scale_fill + shared_theme,
    abs_plot + shared_scale_fill + shared_theme,
    nrow = 2, common.legend = T, legend = "right")
  
  ggsave(
    file = paste(outdir, "sampling_intensity.png", sep = "/"), 
    plot = sampling_intensity_plot)
    
  return(plot_data)
}

