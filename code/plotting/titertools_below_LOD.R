# plot study by country of corresponding author to show data distribution
rm(list = ls())

library(meantiter)
library(titertools)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(grid)
library(gtable)
library(patchwork)
library(ggpubr)

set.seed(100)

# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
figures_dir <- file.path(working_dir, "figures", "over_titer")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")
tables_dir <- file.path(data_dir, "summary_tables")
#----------------------------------------------- set path to save -----------------------------------------------
fileext <- "png"
plot_titers <- TRUE

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))
source(file.path(utility_dir, "plot_since_exposure_functions.R"))

#-------------------------------------  SET TABLE NAME

# here for multiple tables
table_names = list('omicron_folddrops_preprocessed_wSAVE_lessWD.csv')

#for (table_name in table_names){
for (table_name in table_names){
  
  
  # create path to save for each table
  tab_name <- strsplit(table_name, "\\.")[[1]][1]
  
  path_to_save <- file.path(figures_dir, tab_name)
  
  # create subfolder for boxplots
  suppressWarnings(dir.create(path_to_save, recursive = TRUE))
  
  
  table_path <- file.path(google_sheets_dir,table_name)
  
  
  #----------------- load and prepare data
  forest_data <- read.csv(table_path)
  
  
  forest_data <- format_table_forest_plots(forest_data) %>%
    filter(Webplotdigitizer != "y") %>%
    filter(Webplotdigitizer != "y (Titers)")
  
  # select only 2x Vax and BA.1 Omicron. Most data in this group, most relevant occurrence
  forest_data %>% 
    filter(standardise_encounters %in% c("3x Vax", "2x Vax")) %>%
    filter(`Comparator antigen` %in% c("D614G", "WT")) %>%
    filter(OmicronVariant == "BA.1") -> forest_data
  
  
  # set WT titers that are <20 to <2*Value
  # for those also set Omicron values to <2*value
  # for all others set Omicron values to <2*Omicron titer
  # potentially do the same for the ones with ">"
  
  forest_data %>%
    mutate(TitersOmicron_thresh = ifelse(Uncertainty == ">>", ifelse(TitersHAg < 20, "<20", paste0("<", 2*TitersOmicron)), TitersOmicron),
          TitersHAg_thresh = ifelse(Uncertainty ==  ">>" & TitersHAg < 20, "<20", TitersHAg)) %>%
    filter(!is.na(TitersHAg)) -> forest_data_thresh
  
  impute_titers_and_fc <- function(forest_data_thresh){
    
    # now we make titer vector for titertools function
    wt_titers <- forest_data_thresh$TitersHAg_thresh
    ba1_titers <- forest_data_thresh$TitersOmicron_thresh
    
    # titer imputation for wt
    gmt_dist <- titertools::gmt(wt_titers, dilution_stepsize = 0)
    imputed_titers <- titertools::impute_gmt_titers(gmt_dist, wt_titers)
    
    
    ba1_dist <- titertools::gmt(ba1_titers, dilution_stepsize = 0)
    ba1_imputed <- titertools::impute_gmt_titers(ba1_dist, ba1_titers)
    
    forest_data_thresh %>%
      mutate(imputed_wt_titers = as.numeric(imputed_titers),
             imputed_ba1_titers = as.numeric(ba1_imputed),
             imputed_fc = imputed_wt_titers/imputed_ba1_titers,
             log2_imputed_fc = log2(imputed_fc),
             log2_wt = log2(imputed_wt_titers/10),
             log2_ba1 = log2(imputed_ba1_titers/10)) -> forest_data_thresh
    
    return(forest_data_thresh)
    
  }
  
  forest_data_imputed <- rbind(impute_titers_and_fc(forest_data_thresh %>%
                                                      filter(standardise_encounters == "2x Vax")),
                               impute_titers_and_fc(forest_data_thresh %>%
                                                      filter(standardise_encounters == "3x Vax"))) 
  
  
  clean_names <- c("Log2HAg" = "GMT WT",
                  "Log2Omi"  = "GMT BA.1",
                  "log2_ba1" = "GMT BA.1",
                  "log2_imputed_fc" = "FC WT/BA.1",
                  "log2_num_drop" = "FC WT/BA.1",
                  "log2_wt" = "GMT WT")
  
  forest_data_imputed %>%
    mutate(log2_num_drop = log2(`numerical Titre drop`),
           imputed_length = Uncertainty == ">>") %>%
    group_by(standardise_encounters) %>%
    summarize(total_n = length(Sera_details_no_time),
              `Fraction imputed` = round(length(imputed_length[imputed_length])/total_n,2)) -> f
  
  # get mean data of both
  fc_mean <- forest_data_imputed %>%
    mutate(log2_num_drop = log2(`numerical Titre drop`)) %>%
    group_by(standardise_encounters) %>%
    summarize(across(c("log2_num_drop", "Log2HAg", "Log2Omi", "log2_imputed_fc", "log2_wt", "log2_ba1"), Rmisc::CI)) %>%
    mutate(variable = rep(c("upper", "mean", "lower"), 1)) %>%
    ungroup() %>%
    pivot_longer(cols = c("log2_num_drop", "Log2HAg", "Log2Omi", "log2_imputed_fc", "log2_wt", "log2_ba1")) %>%
    mutate(fc_val = name %in% c("log2_imputed_fc", "log2_num_drop")) %>%
    mutate(lin_val = ifelse(fc_val, round(2^value, 2), round(2^value*10, 2))) %>%
    group_by(standardise_encounters, name) %>%
    summarize(clean_val = paste0(lin_val[variable == "mean"], " (", lin_val[variable == "lower"], "; ", lin_val[variable == "upper"], ")")) %>%
    mutate(imputed = name %in% c("log2_imputed_fc", "log2_wt", "log2_ba1"),
           new_names = clean_names[name]) %>%
    select(!name) %>%
    pivot_wider(names_from = "imputed", values_from = "clean_val") %>%
    left_join(., f %>%
                select(!total_n), by = "standardise_encounters")
  
  
  

  
  colnames(fc_mean) <- c("Serum group", "Variable", "LOD/2", "<LOD imputed", "Fraction imputed")
  
  write.csv(fc_mean, file.path(tables_dir, tab_name, "fc_lod2_fc_imputed.csv"))
  
  forest_data_imputed %>%
    filter(Uncertainty == ">>") %>%
    filter(as.numeric(TitersHAg) < 100) %>%
    select(Study, Sourcelink, standardise_encounters, Sera_details_no_time, TitersHAg, imputed_wt_titers, TitersOmicron, imputed_ba1_titers, `numerical Titre drop`, imputed_fc) %>%
    mutate(imputed_fc = round(imputed_fc, 2))-> idvl_study_data
  
  colnames(idvl_study_data) <- c("Study", "Sourcelink", "Serum group", "Serum", "GMT WT", "GMT WT imputed", "GMT BA.1", "GMT BA.1 imputed", "FC WT/BA.1", "FC WT/BA.1 imputed")
  
  # plot imputed and real fc over WT titer
  forest_data_imputed %>%
    select(standardise_encounters, Row_long, Log2HAg, log_fold_change, log2_imputed_fc, Uncertainty) %>%
    mutate(diff_imputed = log_fold_change - log2_imputed_fc) -> fc_long
  
  
  fc_long %>%
    filter(Uncertainty == ">>") %>%
    ggplot(aes(x = Log2HAg, y = diff_imputed, color = standardise_encounters)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    geom_point() + 
    theme_bw() + 
    scale_y_continuous(labels = function(x) ifelse(x < 0, paste0("-", round(2^abs(x), 1), "x"), paste0(round(2^abs(x), 1), "x")),
                       name = "FC LOD/2 / FC imputed",
                       limits = c(-5, 5),
                       breaks = c(-5:5)) + 
    scale_x_continuous(labels = function(x) round(2^x*10, 1),
                       name = "WT GMT") + 
    scale_color_manual(values = c("3x Vax" = "grey60",
                                  "2x Vax" = "grey20"),
                       name = "Serum group") + 
    theme(legend.position = c(0.8, 0.8)) -> fc_plot
  
  
  ggsave(file.path(figures_dir, "fc_imputed_over_wt_gmt.png"), fc_plot, dpi = 300, width = 7, height = 4)
 
}
