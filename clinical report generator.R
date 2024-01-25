library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(readr)
library(readxl)
library(ggpubr)
source('func/report_functions.R')
source('func/calculate_scores.R')

# data input ---------------------
panelInfoPath <- 'data/ISG extended panel - for Nanostring.xlsx'
dataPath <- 'data/nanostring_normalized_counts.csv'
sampleInfoPath <- 'data/sample and physician info.xlsx'

sample_info_all <- left_join(
  readxl::read_excel(sampleInfoPath, sheet = 1),
  readxl::read_excel(sampleInfoPath, sheet = 2),
  by = 'Referring physician'
) %>% mutate(`Referring physician` = if_else(is.na(`Referring physician`), 'unknown', `Referring physician`))

panel_info <- readxl::read_excel(path = panelInfoPath, sheet = 1) %>% filter(type=='Endogenous')
dat <- read_csv(dataPath)

ISG_panel <- na.omit(panel_info$`ISG score`)
NFkb_panel <- na.omit(panel_info$`NF-kB score`)
IFNg_panel <- na.omit(panel_info$`IFN-g Score`)

ISG <- geomean_score(dat, ISG_panel) %>%
  add_column(zscore = (zscore_score(dat, ISG_panel))$zscore)
datFiltered <- dat %>% filter(!batch %in% c('EXP-21-DN5206', 'EXP-21-DN5207','EXP-21-DN5208','EXP-21-DN5209',
                                            'EXP-21-DN5210','EXP-21-DN5211','EXP-21-DN5212','EXP-21-DN5214'))
NFkb <- geomean_score(datFiltered, NFkb_panel) %>%
  add_column(zscore = (zscore_score(datFiltered, NFkb_panel))$zscore)
IFNg <- geomean_score(datFiltered, IFNg_panel) %>%
  add_column(zscore = (zscore_score(datFiltered, IFNg_panel))$zscore)


# generate the pdf report -------------

# patient <- 'ISG-6' # SAVI patient and the Mother in Figure 1B, 1C
# patient <- 'ISG-34' # SAVI patient and the Father in Figure 1E, 1F
patient <- 'ISG-1' # CANDLE patient in Figure 2D
rmarkdown::render('func/ISG_report_template.Rmd',
                  params = list(ISG = ISG,
                                NFkb = NFkb,
                                IFNg = IFNg,
                                panel_info = panel_info,
                                cur_sample = patient,
                                sample_info = sample_info_all %>% filter(`Patient ID` == patient)),
                  # output_file = file.path(getwd(), 'reports/test.pdf'))
                  output_file = file.path(getwd(), '..', 'reports', paste0(patient, '.pdf')))



