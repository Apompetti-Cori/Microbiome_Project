######## RRBS FILTER 75% COVERAGE ########
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  message("loading data")
  #load in data
  rrbs_compile <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_compile"))

  #perform coverage filter only on the gfwt samples, also removing liv_4402 due to bad coverage
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & sample_id != "liv_4402") %>% rownames()
  cov <- getCoverage(rrbs_compile, type = "Cov")
  cov <- cov[,gfwt_samples]
  
  message("gathering loci that pass filter")
  #Keep only loci where there coverage greater than 10 in at least 75% of the samples
  keepLoci <- which(DelayedMatrixStats::rowSums2(cov >= 10) >= ncol(cov)*.75)
  saveRDS(keepLoci, here("results","rds","loci_filter_75.rds"))
  
  message("filter out loci")
  rrbs_compile <- rrbs_compile[keepLoci,]
  
  meth_raw <- getMeth(rrbs_compile, type = "raw")
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("filter for variance")
  #Filter loci for the top 50% most variable
  var <- DelayedMatrixStats::rowVars(meth_raw, na.rm = TRUE)
  o <- order(var, decreasing = TRUE)
  o <- o <= (nrow(meth_raw) * 0.5)
  
  rrbs_compile <- rrbs_compile[o,]
  
  message("saving output")
  HDF5Array::saveHDF5SummarizedExperiment(rrbs_compile, 
                                          dir = here("results","h5","rrbs_compile"), 
                                          replace = FALSE,
                                          verbose = TRUE,
                                          prefix = "filter75")
  
  job::export(c())
}, title = "Filtering methylation matrix 75% coverage")

######## RRBS FILTER 66% COVERAGE ########
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  message("loading data")
  #load in data
  rrbs_compile <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_compile"))

  #perform coverage filter only on the gfwt samples, also removing liv_4402 due to bad coverage
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & sample_id != "liv_4402") %>% rownames()
  cov <- getCoverage(rrbs_compile, type = "Cov")
  cov <- cov[,gfwt_samples]
  
  message("gathering loci that pass filter")
  #Keep only loci where there coverage greater than 10 in at least 66% of the samples
  keepLoci <- which(DelayedMatrixStats::rowSums2(cov >= 10) >= ncol(cov)*.66)
  saveRDS(keepLoci, here("results","rds","loci_filter_66.rds"))
  
  message("filter out loci")
  rrbs_compile <- rrbs_compile[keepLoci,]
  
  meth_raw <- getMeth(rrbs_compile, type = "raw")
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("filter for variance")
  #Filter loci for the top 50% most variable
  var <- DelayedMatrixStats::rowVars(meth_raw, na.rm = TRUE)
  o <- order(var, decreasing = TRUE)
  o <- o <= (nrow(meth_raw) * 0.5)
  
  rrbs_compile <- rrbs_compile[o,]
  
  message("saving output")
  HDF5Array::saveHDF5SummarizedExperiment(rrbs_compile, 
                                          dir = here("results","h5","rrbs_compile"), 
                                          replace = FALSE,
                                          verbose = TRUE,
                                          prefix = "filter66")
  
  job::export(c())
}, title = "Filtering methylation matrix 66% coverage")

######## RRBS FILTER 2 SAMPLE COVERAGE ########
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  message("loading data")
  #load in data
  rrbs_compile <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_compile"))
  cov <- getCoverage(rrbs_compile, type = "Cov")
  
  message("gathering loci that pass filter")
  #Keep only loci where there coverage greater than 10 in at least 75% of the samples
  keepLoci <- which(DelayedMatrixStats::rowSums2(cov >= 10) >= 2)
  saveRDS(keepLoci, here("results","rds","loci_filter_2.rds"))
  
  meth_raw <- getMeth(rrbs_compile, type = "raw")
  seqnames <- GenomeInfoDb::seqnames(rrbs_compile) %>% as.character()
  starts <- BiocGenerics::start(rrbs_compile)
  cpgs <- stringr::str_c(seqnames, starts, starts, sep = ".")
  
  message("filter out loci")
  meth_raw <- meth_raw[keepLoci,]
  
  message("filter for variance")
  #Filter loci for the top 50% most variable
  var <- DelayedMatrixStats::rowVars(meth_raw, na.rm = TRUE)
  o <- order(var, decreasing = TRUE)
  o <- o <= (nrow(meth_raw) * 0.5)
  meth_raw <- meth_raw[o, ]
  
  message("filter samples for too many NA's")
  #filter for too many NA's
  nafilter <- DelayedMatrixStats::colSums2(is.na(meth_raw), useNames = TRUE) >= nrow(meth_raw)*.9
  meth_raw <- meth_raw[,!nafilter]
  
  message("saving output")
  HDF5Array::writeHDF5Array(meth_raw, here("results","h5","rrbs_compile","meth_raw_filter_2.h5"))
  job::export(c())
}, title = "Filtering methylation matrix 2 coverage")
