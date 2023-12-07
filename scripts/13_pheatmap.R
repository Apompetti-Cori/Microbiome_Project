source(here::here("scripts", "01_data_preprocessing.R"))
library(glmnet)
library(caret)
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "smooth")
age_sites <- readRDS(here("results/rds/age_sites.rds"))

#filter for colon samples
se <- se[,se$organ == "COL" & se$microbiome %in% c("gf", "wt")]

#filter for aging sites
se <- subsetByOverlaps(se, age_sites)

#filter for coverage/variance
se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "all", progressbar = TRUE)

meth <- getMeth(se, type = "raw")
meth <- realize_Parallel(meth, cores = 6, nblocks = 3)
meth <- meth[!matrixStats::rowAnyNAs(meth),]

colorcode <- colorRampPalette(c('dodgerblue4', 'snow', 'gold'))(51)
test <- ComplexHeatmap::pheatmap(meth, scale = "none", color = colorcode, 
         annotation_col = top_ha)


