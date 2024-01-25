# geomean ISG score
geomean_score <- function(dat, columnSet){
  #columnSet <- raw$gene_name
  ISG <- dat %>%
    rowwise() %>%
    mutate(geomean = exp(mean(log(c_across(all_of(columnSet)) +1 ))))
  return(ISG)
}


# zscore ISG score
zscore_score <- function(dat, columnSet){
  hc <- dat %>% filter(control) 
  hcmean <- hc %>%
    summarise(across(all_of(columnSet), mean))
  hcstd <- hc %>%
    summarise(across(all_of(columnSet), sd))
  
  zscore <-  dat %>% 
    select(columnSet) %>% 
    apply(.,1, function(x){(x-hcmean)/hcstd}) %>% 
    bind_rows() %>% 
    rowSums()
  ISG <- tibble(Sample = dat$Sample, zscore = zscore)
  return(ISG)
}