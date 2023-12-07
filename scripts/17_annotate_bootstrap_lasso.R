source(here::here("scripts", "01_data_preprocessing.R"))
library(annotatr)
#Select wt microbiome samples and perform cov/var filter
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
se <- se[,se$organ == "SPL" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402"]
se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "all", progressbar = FALSE)

#Filter SE for NA's in meth matrix
meth <- getMeth(se, type = "raw")
meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
se <- se[!matrixStats::rowAnyNAs(meth),]

#Get granges of the cpgs (MUST use the se used in the boostrapping)
granges_0 <- se[cpgs_0,] %>% granges()
granges_90 <- se[cpgs_90,] %>% granges()

annots = c('mm10_cpgs', 'mm10_basicgenes')
annotations <- build_annotations(genome = 'mm10', annotations = annots)

gr_annotated = annotate_regions(
  regions = granges_0,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

