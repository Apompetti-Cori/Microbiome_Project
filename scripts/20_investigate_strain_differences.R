job::empty({
  library(limma)
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Load SE, select SPL, wt, age 4
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402" & se$age == 4]
  granges <- readRDS(here("results/rds/18/granges_spl_strain.rds"))
  se <- IRanges::subsetByOverlaps(se, granges, invert = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Obtain methylation values
  meth <- bsseq::getMeth(se, type = "raw") %>% as.matrix()
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #construct design matrix on strain
  design <- model.matrix(~0 + strain, data = colData(se))
  colnames(design) <- c("strain129svev", "strainLgr5")
  contrast.matrix <- limma::makeContrasts(
    by_strain = strain129svev - strainLgr5,
    levels=design
    )
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- limma::topTable(fit2, number = Inf)
  
  sites <- rownames(results[results$adj.P.Val < .05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/20/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/20/topTable_spl_4.rds"))
  saveRDS(granges, here("results/rds/20/granges_spl_4.rds"))
  
  job::export(c())
}, title = "Strain specific differences in spleen at age 4")

job::empty({
  library(limma)
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Load SE, select SPL, wt, age 4
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("data","se","microbiome_se"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt")]
  se <- filter_rrbs(se, min_coverage = 30, max_coverage = 500, percent_coverage = .75, nblocks = 20, workers = 40, cpg_select = "all", progressbar = FALSE)
  se <- se[,se$microbiome == "wt" & se$age == 4]
  granges <- readRDS(here::here("results/rds/18/granges_spl_max_cov_500.rds"))
  se <- IRanges::subsetByOverlaps(se, granges)
  
  #Get methylation values
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6) 
  rownames(meth) <- mcols(se)$chr_base
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #construct design matrix on strain
  design <- model.matrix(~0 + strain, data = colData(se))
  colnames(design) <- c("strain129svev", "strainc57")
  contrast.matrix <- limma::makeContrasts(
    by_strain = strain129svev - strainc57,
    levels=design
  )
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- limma::topTable(fit2, number = Inf)
  
  sites <- rownames(results[results$adj.P.Val < .05,])
  granges <- subset(se, chr_base %in% sites) %>% granges()
  
  #Save results
  dir.create(here("results/rds/20/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/20/topTable_spl_4_max_cov_500.rds"))
  saveRDS(granges, here("results/rds/20/granges_spl_4_max_cov_500.rds"))
  
  job::export(c())
}, title = "Strain specific differences in spleen at age 4 max cov 500")


job::empty({
  library(limma)
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Load SE, select SPL, wt, age 4
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "smooth")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402" & se$age == 4]
  granges <- readRDS(here("results/rds/18/granges_spl_smooth.rds"))
  se <- IRanges::subsetByOverlaps(se, granges, invert = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "smooth")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Obtain methylation values
  meth <- bsseq::getMeth(se, type = "smooth") %>% as.matrix()
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #construct design matrix on strain
  design <- model.matrix(~0 + strain, data = colData(se))
  colnames(design) <- c("strain129svev", "strainc57")
  contrast.matrix <- limma::makeContrasts(
    by_strain = strain129svev - strainc57,
    levels=design
  )
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- limma::topTable(fit2, number = Inf)
  
  sites <- rownames(results[results$adj.P.Val < .05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/20/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/20/topTable_spl_4_smooth.rds"))
  saveRDS(granges, here("results/rds/20/granges_spl_4_smooth.rds"))
  
  job::export(c())
}, title = "Strain specific differences in spleen at age 4 (using smoothed meth)")

job::empty({
  library(limma)
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Load SE, select SPL, wt, age 4
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("data","se","microbiome_se"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402" & se$age == 4]
  granges <- readRDS(here("results/rds/18/granges_spl_max_cov_Inf.rds"))
  se <- IRanges::subsetByOverlaps(se, granges, invert = FALSE)
  se <- annotate_rrbs(se)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #construct design matrix on strain
  design <- model.matrix(~0 + strain, data = colData(se))
  colnames(design) <- c("strain129svev", "strainc57")
  contrast.matrix <- limma::makeContrasts(
    by_strain = strain129svev - strainc57,
    levels=design
  )
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- limma::topTable(fit2, number = Inf)
  
  sites <- rownames(results[results$adj.P.Val < .05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/20/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/20/topTable_spl_4_max_cov_Inf.rds"))
  saveRDS(granges, here("results/rds/20/granges_spl_4_max_cov_Inf.rds"))
  
  job::export(c())
}, title = "Strain specific differences in spleen at age 4 (using max cov Inf)")

job::empty({
  library(limma)
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Load SE, select SPL, wt, age 4
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("data","se","microbiome_se"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402" & se$age == 4]
  granges <- readRDS(here("results/rds/18/granges_spl_max_cov_Inf_cgi.rds"))
  se <- IRanges::subsetByOverlaps(se, granges, invert = FALSE)
  se <- annotate_rrbs(se)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  rownames(meth) <- mcols(se)$chr_base 
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #construct design matrix on strain
  design <- model.matrix(~0 + strain, data = colData(se))
  colnames(design) <- c("strain129svev", "strainc57")
  contrast.matrix <- limma::makeContrasts(
    by_strain = strain129svev - strainc57,
    levels=design
  )
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- limma::topTable(fit2, number = Inf)
  
  sites <- rownames(results[results$adj.P.Val < .05,])
  granges <- subset(se, chr_base %in% sites) %>% granges()
  
  #Save results
  dir.create(here("results/rds/20/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/20/topTable_spl_4_max_cov_Inf_cgi.rds"))
  saveRDS(granges, here("results/rds/20/granges_spl_4_max_cov_Inf_cgi.rds"))
  
  job::export(c())
}, title = "Strain specific differences in spleen at age 4 (using max cov Inf cgi only)")

job::empty({
  library(limma)
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Load SE, select SPL, wt, age 4
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("data","se","microbiome_se"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402" & se$age == 4]
  granges <- readRDS(here("results/rds/18/granges_spl_max_cov_500.rds"))
  se <- IRanges::subsetByOverlaps(se, granges, invert = FALSE)
  se <- annotate_rrbs(se)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #construct design matrix on strain
  design <- model.matrix(~0 + strain, data = colData(se))
  colnames(design) <- c("strain129svev", "strainc57")
  contrast.matrix <- limma::makeContrasts(
    by_strain = strain129svev - strainc57,
    levels=design
  )
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- limma::topTable(fit2, number = Inf)
  
  sites <- rownames(results[results$adj.P.Val < .05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/20/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/20/topTable_spl_4_max_cov_500.rds"))
  saveRDS(granges, here("results/rds/20/granges_spl_4_max_cov_500.rds"))
  
  job::export(c())
}, title = "Strain specific differences in spleen at age 4 (using max cov 500)")

job::empty({
  library(limma)
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Load SE, select COL, wt
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "COL" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "cgi", progressbar = FALSE)
  
  
  meth <- getMeth(se, type = "raw") %>% as.matrix()
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Obtain methylation values
  meth <- bsseq::getMeth(se[,se$age == 4], type = "raw") %>% as.matrix()
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #construct design matrix on strain
  design <- model.matrix(~0 + strain, data = colData(se))
  colnames(design) <- c("strain129svev", "strainLgr5")
  contrast.matrix <- limma::makeContrasts(
    by_strain = strain129svev - strainLgr5,
    levels=design
  )
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- limma::topTable(fit2, number = Inf)
  
  sites <- rownames(results[results$adj.P.Val < .05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/20/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/20/topTable_col_4.rds"))
  saveRDS(granges, here("results/rds/20/granges_col_4.rds"))
  
  job::export(c())
}, title = "Strain specific differences in colon at age 4")


job::empty({
  library(limma)
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Load SE, select LIV, wt
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "LIV" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402"]
  granges <- readRDS(here("results/rds/18/granges_liv_strain.rds"))
  se <- IRanges::subsetByOverlaps(se, granges, invert = FALSE)
  
  meth <- getMeth(se, type = "raw") %>% as.matrix()
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Obtain methylation values
  meth <- bsseq::getMeth(se[,se$age == 4], type = "raw") %>% as.matrix()
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #construct design matrix on strain
  design <- model.matrix(~0 + strain, data = colData(se))
  colnames(design) <- c("strain129svev", "strainLgr5")
  contrast.matrix <- limma::makeContrasts(
    by_strain = strain129svev - strainLgr5,
    levels=design
  )
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  results <- limma::topTable(fit2, number = Inf)
  
  sites <- rownames(results[results$adj.P.Val < .05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/20/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/20/topTable_liv_4.rds"))
  saveRDS(granges, here("results/rds/20/granges_liv_4.rds"))
  
  job::export(c())
}, title = "Strain specific differences in liver at age 4")