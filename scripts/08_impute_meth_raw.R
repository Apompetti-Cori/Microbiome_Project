job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  set.seed(123)
  
  #Configure delayed array options
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  DelayedArray::setAutoRealizationBackend(BACKEND="HDF5Array")
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #summarize missingness before imputation
  message("nrow: ", nrow(meth_raw))
  message("complete cases before: ", nrow(meth_raw) - sum(DelayedMatrixStats::rowAnyNAs(meth_raw)))
  
  #impute row means through block processing and rbind after
  message("imputing row means")
  imputes <- blockApply(meth_raw, grid = rowAutoGrid(meth_raw, nrow = 10000),
                       function(block) {block.impute.row_mean(block)}, 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 5, progressbar = TRUE))
  meth_raw_impute_rm <- Reduce(rbind, imputes)
  
  #summarize missingness after imputation (should be same as nrow)
  #message("complete cases after:", nrow(meth_raw_impute_rm) - sum(DelayedMatrixStats::rowAnyNAs(meth_raw_impute_rm)))
  
  #save imputed matrix
  message("saving output")
  HDF5Array::writeHDF5Array(meth_raw_impute_rm, 
                            here("results","h5","rrbs_combine","meth_raw_filter_75_impute_rm.h5"), 
                            with.dimnames = TRUE, name = "meth", verbose = TRUE)
  
  DelayedArray:::set_verbose_block_processing(verbose = FALSE)
  job::export(c())
}, title = "Row means imputation")
