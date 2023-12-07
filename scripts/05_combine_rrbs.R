job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  DelayedArray::setAutoRealizationBackend(BACKEND="HDF5Array")
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  message("loading data")
  rrbs_org <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_org"), prefix = "mm10")
  rrbs_gfwt <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  
  message("combining bsseq objects")
  rrbs_compile <- bsseq::combineList(rrbs_gfwt, rrbs_org)
  
  message("saving output")
  HDF5Array::saveHDF5SummarizedExperiment(rrbs_compile, 
                                          dir = here("results","h5","rrbs_compile"), 
                                          replace = TRUE,
                                          verbose = TRUE)
  job::export(c())
}, title = "Combine RRBS")


job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  DelayedArray::setAutoRealizationBackend(BACKEND="HDF5Array")
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  message("loading data")
  rrbs_new <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("data/se/0827"))
  rrbs_old <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("data/se/rrbs_compile"))
  
  message("combining bsseq objects")
  rrbs_compile <- bsseq::combineList(rrbs_old, rrbs_new)
  
  message("saving output")
  HDF5Array::saveHDF5SummarizedExperiment(rrbs_compile, 
                                          dir = here("results","h5","rrbs_compile"), 
                                          replace = FALSE,
                                          prefix = "v2",
                                          verbose = TRUE)
  job::export(c())
}, title = "Combine RRBS")