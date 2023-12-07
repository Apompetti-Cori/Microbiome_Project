########### CONDUCT PCA WITH COMPLETE CASES ############

########### CONDUCT PCA WITH ALL SAMPLES MINUS LIV_4402 ############
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in se
  rrbs_compile <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_compile"), prefix = "filter75")
  meth_raw <- getMeth(rrbs_compile, type = "raw")
  
  #keep everything but liv_4402
  samples <- md_rrbs %>% filter(sample_id != "liv_4402") %>% rownames()
  meth_raw <- meth_raw[,samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep only complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  message("complete cases:",  nrow(meth_raw))
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save ouput
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_all_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA 75")

########### CONDUCT PCA WITH GFWT SAMPLES (COLON) ############
#### BOTH STRAINS ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt LIV samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & sample_id != "liv_4402" & organ %in% c("DCOL","PCOL","COL")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep only complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_colon_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA colon 75")

#### C57BL6 ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt LIV samples within BL6
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" & 
                                                 organ %in% c("DCOL","PCOL","COL") & 
                                                 strain %in% c("Lgr5?EGFP?IRES?CreERT2",
                                                               "c57bl6")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol all: ", ncol(meth_raw))
  
  #keep only complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  message("ncol complete: ", ncol(meth_raw))
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_colon_bl6_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA colon 75 BL6")

#### 129SVEV ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt LIV samples within 129svev
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" & 
                                                 organ %in% c("DCOL","PCOL","COL") & 
                                                 strain %in% c("129svev")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol all: ", ncol(meth_raw))
  
  #keep only complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  message("ncol complete: ", ncol(meth_raw))
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_colon_129_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA colon 75 129svev")

########### CONDUCT PCA WITH GFWT SAMPLES (LIVER MINUS LIV_4402) ############
#### BOTH STRAINS ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt LIV samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" & 
                                                 organ == "LIV") %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_liver_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA liver 75")

#### C57BL6 ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt LIV samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" &
                                                 organ == "LIV" & 
                                                 strain %in% c("Lgr5?EGFP?IRES?CreERT2",
                                                               "c57bl6")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_liver_bl6_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA liver 75 BL6")

#### 129SVEV ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt LIV samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" &
                                                 organ == "LIV" & 
                                                 strain %in% c("129svev")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_liver_129_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA liver 75 129")

########### CONDUCT PCA WITH GFWT SAMPLES (SPLEEN) ############
#### BOTH STRAINS ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt SPL samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & sample_id != "liv_4402" & organ == "SPL") %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca 
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_spleen_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA spleen 75")

#### C57BL6 ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt SPL samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" & 
                                                 organ == "SPL" & 
                                                 strain %in% c("Lgr5?EGFP?IRES?CreERT2",
                                                               "c57bl6")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca 
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_spleen_bl6_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA spleen 75 BL6")

#### 129SVEV ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt SPL samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" & 
                                                 organ == "SPL" & 
                                                 strain %in% c("129svev")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca 
  message("conducting pca")
  pca <- PCAtools::pca(meth_raw, metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_spleen_129_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA spleen 75 129")

########### CONDUCT PCA WITH GFWT SAMPLES (SMALL INTESTINE) ############
#### BOTH STRAINS
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt SI samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" & 
                                                 organ %in% c("LSI", "USI")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca 
  message("conducting pca")
  pca <- PCAtools::pca(na.omit(meth_raw), metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_si_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA SI 75")

#### C57BL6 ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt SI samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" & 
                                                 organ %in% c("LSI", "USI") & 
                                                 strain %in% c("Lgr5?EGFP?IRES?CreERT2",
                                                               "c57bl6")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca 
  message("conducting pca")
  pca <- PCAtools::pca(na.omit(meth_raw), metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_si_bl6_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA SI 75 BL6")

#### 129SVEV THESE DON'T EXIST IN SI YET MAYBE ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep only gfwt SI samples
  gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & 
                                                 sample_id != "liv_4402" & 
                                                 organ %in% c("LSI", "USI") & 
                                                 strain %in% c("129svev")) %>% rownames()
  meth_raw <- meth_raw[,gfwt_samples]
  
  message("ncol: ", ncol(meth_raw))
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca 
  message("conducting pca")
  pca <- PCAtools::pca(na.omit(meth_raw), metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","PCA","pca_75_si_129_cc.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA SI 75 129")

########### CONDUCT PCA WITH IMPUTED VALUES ############
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75_impute_rm.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(na.omit(meth_raw), metadata = md_rrbs[colnames(meth_raw),])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","pca_75_all_rm.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA")

job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #load in h5 file
  message("loading data")
  h5file <- here("results","h5","rrbs_combine","meth_raw_filter_75_impute_rm.h5")
  meth_raw <- HDF5Array::HDF5Array(filepath = h5file, name = "meth")
  
  #keep complete cases
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  
  #keep samples belonging to colon
  col_samples <- md_rrbs %>% filter(organ == c("COL", "DCOL", "PCOL") & batch != "Organoid") %>% rownames()
  meth_raw_col_indx <- colnames(meth_raw)[colnames(meth_raw) %in% col_samples]
  
  #conduct pca
  message("conducting pca")
  pca <- PCAtools::pca(na.omit(meth_raw[,meth_raw_col_indx]), metadata = md_rrbs[meth_raw_col_indx,])
  
  #save output
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","pca_75_colon_rm.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA colon")

job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  message("loading data")
  meth_raw <- trqwe::mcreadRDS(here::here("results", "rds", "meth_raw_filter_75_rowmeans.rds"), mc.cores = 4)
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  message("conducting pca")
  liv_samples <- md_rrbs %>% filter(organ == c("LSI", "USI")) %>% rownames()
  meth_raw_liv_indx <- colnames(meth_raw)[colnames(meth_raw) %in% liv_samples]
  pca <- PCAtools::pca(na.omit(meth_raw[, meth_raw_liv_indx]), metadata = md_rrbs[meth_raw_liv_indx,])
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","pca_75_liver_rm.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA liver")

job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  message("loading data")
  meth_raw <- trqwe::mcreadRDS(here::here("results", "rds", "meth_raw_filter_75_rowmeans.rds"), mc.cores = 4)
  meth_raw <- meth_raw[!DelayedMatrixStats::rowAnyNAs(meth_raw),]
  message("conducting pca")
  spl_samples <- md_rrbs %>% filter(organ == "SPL") %>% rownames()
  meth_raw_spl_indx <- colnames(meth_raw)[colnames(meth_raw) %in% spl_samples]
  pca <- PCAtools::pca(na.omit(meth_raw[, meth_raw_spl_indx]), metadata = md_rrbs[meth_raw_spl_indx,])
  message("saving output")
  trqwe::mcsaveRDS(pca, here("results","rds","pca_75_spleen_rm.rds"), mc.cores = 4)
  job::export(c())
}, title = "PCA spleen")