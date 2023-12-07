source(here::here("scripts", "01_data_preprocessing.R"))
library(glmnet)
library(caret)
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
age_sites <- readRDS(here("results/rds/age_sites.rds"))

#filter for colon samples
se <- se[,se$organ == "COL" & se$microbiome %in% c("gf", "wt")]

#filter for aging sites
se <- subsetByOverlaps(se, age_sites)

#filter for coverage/variance
se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "cgi", progressbar = TRUE)

#get methylation values
meth <- getMeth(se, type = "raw")
meth <- realize_Parallel(meth, cores = 6, nblocks = 3)

