#### Correlations across all sites in spleen ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, workers = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get just the wt samples
  se <- se[,se$microbiome %in% c("wt")]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~age + strain, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  results <- limma::topTable(fit, coef = 2, number = Inf)
  
  #Convert sites to granges
  sites <- rownames(results[results$adj.P.Val < 0.1,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_spl.rds"))
  saveRDS(granges, here("results/rds/18/granges_spl.rds"))
  
  job::export(c())
}, title = "linear model across all sites in Spleen")

#### Correlations across all sites in spleen using min cov 30 max cov 500 75%####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("data","se","microbiome_se"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, min_coverage = 30, max_coverage = 500, percent_coverage = .75, nblocks = 20, workers = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter meth for just wt samples
  se <- se[,se$microbiome == "wt"]
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  rownames(meth) <- mcols(se)$chr_base
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~age + strain, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  results <- limma::topTable(fit, coef = 2, number = Inf)
  
  #Convert sites to granges
  sites <- rownames(results[results$adj.P.Val < 0.1,])
  granges <- subset(se, chr_base %in% sites) %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_spl_max_cov_500.rds"))
  saveRDS(granges, here("results/rds/18/granges_spl_max_cov_500.rds"))
  
  job::export(c())
}, title = "linear model across all sites in Spleen max cov 500")

#### Correlations across all sites in spleen using min cov 30 max cov Inf 75%####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, min_coverage = 30, max_coverage = Inf, percent_coverage = .75, nblocks = 20, workers = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter for wt samples
  se <- se[,se$microbiome == "wt"]
  meth <- getMeth(se, type = "raw")
  meth <- meth %>% as.matrix()
  rownames(meth) <- mcols(se)$chr_base
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~ 0 + age + strain, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  results <- limma::topTable(fit, coef = 1, number = Inf)

  #Convert sites to granges
  sites <- rownames(results[results$adj.P.Val < 0.1,])
  granges <- subset(se, chr_base %in% sites) %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_spl_max_cov_Inf.rds"))
  saveRDS(granges, here("results/rds/18/granges_spl_max_cov_Inf.rds"))
  
  job::export(c())
}, title = "linear model across all sites in Spleen max cov Inf")

#### Correlations across all sites in spleen using min cov 30 max cov Inf 75% CpG Islands only####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, min_coverage = 30, max_coverage = Inf, percent_coverage = .75, nblocks = 20, workers = 40, cpg_select = "cgi", progressbar = FALSE)
  
  #Filter for wt samples
  se <- se[,se$microbiome == "wt"]
  meth <- getMeth(se, type = "raw")
  meth <- meth %>% as.matrix()
  rownames(meth) <- mcols(se)$chr_base
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~ 0 + age + strain, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  results <- limma::topTable(fit, coef = 1, number = Inf)
  
  #Convert sites to granges
  sites <- rownames(results[results$adj.P.Val < 0.1,])
  granges <- subset(se, chr_base %in% sites) %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_spl_max_cov_Inf_cgi.rds"))
  saveRDS(granges, here("results/rds/18/granges_spl_max_cov_Inf_cgi.rds"))
  
  job::export(c())
}, title = "linear model across all sites in Spleen max cov Inf CpG Islands only")

#### Correlations across all sites in spleen (using smoothed values) ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "smooth")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, workers = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "smooth")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get just the wt samples
  se <- se[,se$microbiome %in% c("wt")]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "smooth")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~age + strain, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  results <- limma::topTable(fit, coef = 2, number = Inf)
  
  #Convert sites to granges
  sites <- rownames(results[results$adj.P.Val < 0.05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_spl_smooth.rds"))
  saveRDS(granges, here("results/rds/18/granges_spl_smooth.rds"))
  
  job::export(c())
}, title = "linear model across all sites in Spleen (using smoothed meth vals)")

#### Correlations across all sites in spleen removing 129svev ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, workers = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get just the wt samples
  se <- se[,se$microbiome %in% c("wt")]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~age + strain, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  results <- limma::topTable(fit, coef = 2, number = Inf)
  
  #Convert sites to granges
  sites <- rownames(results[results$adj.P.Val < 0.05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_spl_bl6.rds"))
  saveRDS(granges, here("results/rds/18/granges_spl_bl6.rds"))
  
  job::export(c())
}, title = "linear model across all sites in Spleen c57bl6 only")

#### Correlations across cpgi sites in spleen ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, workers = 40, cpg_select = "cgi", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get just the wt samples
  se <- se[,se$microbiome %in% c("wt")]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~age + strain, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  #Convert sites to granges
  results <- limma::topTable(fit, coef = 2, number = Inf)
  sites <- rownames(results[results$adj.P.Val < 0.05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_spl_cgi.rds"))
  saveRDS(granges, here("results/rds/18/granges_spl_cgi.rds"))
  
  job::export(c())
}, title = "linear model across cgi sites in Spleen")

#### Correlations across all sites in liver ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "LIV" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, workers = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get just the wt samples
  se <- se[,se$microbiome %in% c("wt")]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~ 0 + age + strain, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  results <- limma::topTable(fit, coef = 1, number = Inf)
  
  #Convert sites to granges
  sites <- rownames(results[results$adj.P.Val < 0.1,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_liv_strain.rds"))
  saveRDS(granges, here("results/rds/18/granges_liv_strain.rds"))
  
  job::export(c())
}, title = "linear model across all sites in Liver controlling for strain")


#### Correlations across all sites in orginal samples ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #Select gf/wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  #obtain the samples Himani used for aging sites
  samples <- read.csv(here("data/metadata/annotated_perm_intestine44_age.tsv"), sep = "\t") %>% 
    colnames() %>% 
    grep(pattern = ".cov.gz", value = TRUE) %>% 
    tools::file_path_sans_ext() %>% 
    tools::file_path_sans_ext()
  
  #filter se for Himani's samples
  se <- se[,samples]
  
  #filter for cov/variance
  se <- filter_rrbs(se, nblocks = 20, workers = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, workers = 12, nblocks = 6)
  
  #convert beta values to m values
  mvals <- beta2m(meth)
  
  #create design matrix
  design <- model.matrix(~age, data = colData(se))
  
  #fit mvals on design matrix
  fit <- limma::lmFit(mvals, design)
  
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  
  #Convert sites to granges
  results <- limma::topTable(fit, coef = 2, number = Inf)
  sites <- rownames(results[results$adj.P.Val < 0.05,]) %>% as.numeric()
  granges <- se[sites,] %>% granges()
  
  #Save results
  dir.create(here("results/rds/18/"), showWarnings = FALSE)
  saveRDS(results, here("results/rds/18/topTable_himani.rds"))
  saveRDS(granges, here("results/rds/18/granges_himani.rds"))
  
  job::export(c())
}, title = "linear model across cgi sites in Orginal Samples")
