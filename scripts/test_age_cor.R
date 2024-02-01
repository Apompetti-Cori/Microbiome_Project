library(tidyverse)
se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs/"), prefix = "isee_filter_side")

dm_res <- contrastResults(se, "cons_gf") %>% as.data.frame()
rownames(dm_res) <- dm_res$chr_base
