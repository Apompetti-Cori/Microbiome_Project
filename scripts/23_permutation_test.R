se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("data","se","rrbs_compile"))




meth <- getMeth(se, type = "raw") 
