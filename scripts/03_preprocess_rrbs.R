# PROCESS RRBS ORG
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  
  covfiles <- paste(here("data/covfiles/"), md_rrbs %>% filter(batch == "Organoid") %>% pull("filename"), sep = "")
  colData <- md_rrbs %>% filter(batch == "Organoid")
  rrbs_org <- bsseq::read.bismark(files = covfiles,
                              verbose = TRUE,
                              colData = colData,
                              BACKEND = "HDF5Array",
                              dir = here("results", "h5", "rrbs_org"),
                              replace = TRUE,
                              BPPARAM = BiocParallel::MulticoreParam(workers = 16, progressbar = TRUE)
  )
  
  #Convert RRBS Org to mm10
  
  library(rtracklayer)
  mm9tomm10 <- import.chain(here("data", "chainfiles","mm9ToMm10.over.chain"))
  
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
}, title = "Process RRBS ORG")

# PROCESS RRBS GFWT
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  
  md_rrbs <- read.csv(here(""))
  
  covfiles <- paste(here("data/covfiles/"), md_rrbs %>% filter(sampletype != "Organoid") %>% pull(filename), sep = "")
  colData <- md_rrbs %>% filter(sampletype != "Organoid")
  rrbs_gfwt <- bsseq::read.bismark(files = covfiles,
                                   colData = colData,
                                   verbose = TRUE,
                                   BACKEND = "HDF5Array",
                                   dir = here("data/se/0827"),
                                   replace = TRUE,
                                   BPPARAM = BiocParallel::MulticoreParam(workers = 60, progressbar = TRUE)
  )
  
  message("converting seqnames")
  rrbs_gfwt <- renameSeqlevels(rrbs_gfwt, c(MT="M"))
  rrbs_gfwt <- renameSeqlevels(rrbs_gfwt, paste("chr",seqlevels(rrbs_gfwt),sep=""))
  
  message("saving changes")
  HDF5Array::saveHDF5SummarizedExperiment(x = rrbs_gfwt, 
                                          dir = here("results","h5","rrbs_gfwt"), 
                                          replace = FALSE,
                                          prefix = "update")
  
  job::export(c())
}, title = "Process RRBS GFWT")