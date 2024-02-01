message("preprocessing")
HDF5Array::setHDF5DumpDir(here::here("results","h5","HDF5DumpDir"))
DelayedArray::setAutoRealizationBackend(BACKEND="HDF5Array")
suppressPackageStartupMessages({
  library(bsseq)
  library(tidyverse)
  library(PCAtools) 
  library(ggplot2) 
  library(plotly) 
  library(here) 
  library(ggeasy) 
  library(rhdf5) 
  library(HDF5Array) 
  library(data.table)
  library(DelayedArray)
  library(DelayedMatrixStats) 
})

if("functions" %in% search()){
  detach("functions")
}

functions <- new.env()
sys.source(here::here("scripts/functions/functions.R"), envir = functions)
attach(functions)
rm(functions)