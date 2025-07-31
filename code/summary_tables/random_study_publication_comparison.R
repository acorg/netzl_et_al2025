rm(list = ls())


library(tidyverse)
library(readxl)

# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
figures_dir <- file.path(working_dir, "figures", "by_country")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "plot_over_time_functions.R"))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))

set.seed(100)
#-------------------------------------  Read in studies table
studies <- read.csv(file.path(google_sheets_dir, "Netzl et al. - Collected Omicron antigenic data.csv")) %>%
  select(Study, Sourcelink, Country, year) %>%
  unique()


selected_studies <- studies[sample(1:85, 8, replace = FALSE),]

publication_links <- c("https://doi.org/10.1016/j.chom.2022.07.002", #new published study
                       "https://doi.org/10.1080/22221751.2022.2099305", #new published study
                       "https://www.cell.com/cell/pdf/S0092-8674(21)01578-6.pdf?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421015786%3Fshowall%3Dtrue", # alredy published paper in manuscript, no new data
                       "https://doi.org/10.3389/fimmu.2022.946318", #new published study
                       "https://doi.org/10.1038/s41421-022-00375-5", #new published study
                       "https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(22)00060-X/fulltext", #already published, no new data
                       "https://doi.org/10.1016/j.chom.2022.09.018", #new published study
                       "https://doi.org/10.1056/nejmc2119426" #new published study
)


# note on stiasny: Additional samples resulted in change, others did not change at all
# last one had different group labels and some additional data

# read in preprint and published data
preprint_data <- read_xlsx(file.path(google_sheets_dir, "omicron_folddrops_preprocessed_wSAVE_lessWD_publication_comparison.xlsx"), sheet = 1) %>%
  mutate(Data = "Preprint") 

published_data <- read_xlsx(file.path(google_sheets_dir, "omicron_folddrops_preprocessed_wSAVE_lessWD_publication_comparison.xlsx"), sheet = 2) %>%
  mutate(Data = "Publication") 

comb_data <- rbind(preprint_data, 
                   published_data)


# make WT titers to long
comb_data <- rbind(comb_data,
                   comb_data %>%
                     mutate(OmicronVariant = `Comparator antigen`,
                            TitersOmicron = TitersHAg)) %>%
  select(!TitersHAg) %>%
  select(!`Comparator antigen`) %>%
  unique()
  
# make table that shows both titers and number of sera
comb_data %>%
  select(Study, Sourcelink, `Sera details long`, `Number of sera`, OmicronVariant, TitersOmicron, Data) %>%
  unique() %>%
  group_by(Study) %>%
  mutate(Sources = paste(unique(c(Sourcelink[Data == "Preprint"], Sourcelink[Data == "Publication"])), collapse = " \n ")) %>%
  select(Study, Sources, `Sera details long`, `Number of sera`, OmicronVariant, TitersOmicron, Data) %>%
  pivot_wider(names_from = c("Data"), values_from = c("Number of sera", "TitersOmicron")) %>%
  relocate(Study, Sources, `Sera details long`, `Number of sera_Preprint`,
           `Number of sera_Publication`,
           OmicronVariant,
           TitersOmicron_Preprint,
           TitersOmicron_Publication) %>%
  mutate(OmicronVariant = factor(OmicronVariant, 
                                    levels = c("WT", "D614G", "BA.1", "BA.1.1", "BA.2", "BA.2.12.1",
                                               "BA.3", "BA.4/5", "BA.2.38", "BA.2.75", "BA.2.76"))) %>%
  arrange(Study, `Sera details long`, OmicronVariant)-> comb_wide


# substitute repeat values with empty line for table printing
duplicated_lines <- duplicated(comb_wide$`Sera details long`) & duplicated(comb_wide$Sources)

comb_wide[duplicated_lines, c("Study", "Sources", "Sera details long", "Number of sera_Preprint", "Number of sera_Publication")] <- NA


comb_wide[duplicated(comb_wide$Sources), c("Study", "Sources")] <- NA

write.csv(comb_wide, file.path(google_sheets_dir, "preprint_publication_subset_comparison.csv"))
