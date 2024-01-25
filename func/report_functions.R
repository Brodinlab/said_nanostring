library(dplyr)
library(purrr)
library(kableExtra)
library(ggplot2)
library(ggrepel)


neat_print_table <- function(df,col_per_row, col_names=NA, caption = NULL, ...){
  # the extra arguments will be passed to kable_styling
  split(1:ncol(df), sort(rep_len(1:ceiling(ncol(df)/col_per_row), ncol(df)))) %>%
    map(~dplyr::select(df, .)) %>%
    map(kbl, booktabs = T, col.names = col_names, caption = caption) %>%
    map(kable_styling, latex_options = c('hold_position', 'striped'), ...) %>%
    walk(print)
}

prepare_score_table <- function(df, panel){
  tab <- rbind(
    t(df[,c('geomean', 'zscore')]),
    t(df[,panel])
  )
  colnames(tab) <- make.unique(paste0(df$`Subject ID`, '_Visit', df$Visit))
  tab <- round(tab)
  if (dim(tab)[2] > 1){
    tab <- tab[,gtools::mixedsort(colnames(tab))]
  } # sort column when there's more than 1 column. Otherwise the column name is lost 
  tab <- as.data.frame(tab)
  return(tab)
}

score_scatter_plot <- function(df){
  set.seed(42)
  ggplot(df, aes(x = Group, y = geomean, color=Group)) +
    geom_jitter(width = 0.25) +
    geom_text_repel(aes(label = label),
                    force = 2,
                    max.overlaps = 20,
                    segment.color = 'transparent') +
    theme_bw()  +
    xlab('') +
    theme(legend.position='none',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_y_continuous(trans='log1p')
}

prepare_timeline_table <- function(df, current_sample, sample_info){
  d <- df %>% filter(`Subject ID` == current_sample) %>%
    mutate(label = paste(`Subject ID`, Visit, sep='_')) %>%
    left_join(sample_info %>% 
                mutate(label = paste(`Patient ID`, Visit, sep='_')) %>%
                select(all_of(c('Date of sampling', 'label'))), 
              by='label') %>%
    mutate(`Date of sampling` = as.Date(`Date of sampling`),
           connection_group = if_else(Group == 'Patients' & !grepl('M|F', Visit), TRUE, FALSE)) %>%
    filter(!grepl('stim', Visit)) %>%
    arrange(`Date of sampling`)
  return(d)
}

score_timeline_plot <- function(df, healthy_range){
  set.seed(42)
  x_min = min(df$`Date of sampling`, na.rm = T) -100
  x_max = max(df$`Date of sampling`, na.rm = T) +100
  ggplot(df, aes(x=`Date of sampling`, y=geomean)) + 
    geom_rect(aes(xmin = x_min, xmax = x_max, ymin = healthy_range[1], ymax = healthy_range[2]),
              fill = "lightblue", alpha=0.05) +
    geom_point() + 
    geom_path(data = . %>% filter(connection_group), aes(group=connection_group)) +
    scale_x_date(date_labels = "%Y %b %d", breaks = unique(df$`Date of sampling`)) + 
    geom_text_repel(aes(label=label), size = 2.5) +
    # scale_y_continuous(trans='log1p') + 
    annotate(geom="text", x=x_min, y=healthy_range[2], hjust = 0, vjust = 1, label="Range of healthy controls",
             color="blue")+
    theme_bw() +
    theme(aspect.ratio=3/8, 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.minor.x = element_blank(),
          panel.background = element_blank())
}








