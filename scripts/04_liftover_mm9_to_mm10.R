job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  library(rtracklayer)
  mm9tomm10 <- import.chain(here("data", "chainfiles","mm9ToMm10.over.chain"))
  message("loading data")
  rrbs_org <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_org"))
  
  mcols(rrbs_org)$index <- 1:length(rrbs_org)
  mm9 <- rowRanges(rrbs_org)
  mm9$index <- 1:length(mm9)
  
  #liftover mm9 to mm10
  mm10 <- mm9 %>% liftOver(mm9tomm10) %>% unlist()
  
  # Subset based on regions that lifted over
  rrbs_org.mm10 <- rrbs_org[which(mm10$index %in% mm9$index)]
  
  rowRanges(rrbs_org.mm10) <- mm10
  
  print(glue::glue("{length(rrbs_org.mm10)} out of {length(rrbs_org)} were lifted over"))
  print(glue::glue("{length(rrbs_org) - length(rrbs_org.mm10)} did not liftOver"))
  
  genome(rrbs_org.mm10) <- "mm10"
  
  HDF5Array::saveHDF5SummarizedExperiment(x = rrbs_org.mm10, 
                                          dir = here("results","h5","rrbs_org"), 
                                          replace = TRUE,
                                          prefix = "mm10")
  job::export(c())
}, title = "Liftover mm9 to mm10")
