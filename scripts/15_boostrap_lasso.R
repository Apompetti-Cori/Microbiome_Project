#### Bootstrap select spleen aging loci ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  library(glmnet)
  library(boot)
  library(data.table)
  library(doMC)
  
  #Select wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
  meth <- t(meth)
  
  #Get actual age of samples
  age <- colData(se)[,"age"]
  
  #Bind age to meth values
  meth_age <- cbind(age,meth)
  
  #Define function that will be run in bootstrapping
  lasso_boot <- function(){
    
    rsamples <- as.data.table(colData(se))[, .SD[sample(.N, replace = TRUE)], by = age][, sample_id]
    
    rmeth_age <- meth[rsamples,]
    
    fit <- cv.glmnet(x = rmeth_age[,-1], y = rmeth_age[,1], alpha = 1)
    
    return(which(coef(fit) > 0))
  }
  
  #Define how many boots and cores will be used for bootstrapping
  boots <- 1000
  cores <- 8
  
  boot_results <- parallel::mclapply(
    1:boots, 
    function(i) suppressMessages(lasso_boot()),
    mc.cores = cores
  )
  
  job::export(c(boot_results))
},title = "Bootstrap select spleen aging loci")

#### Bootstrap select SI aging loci ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  library(glmnet)
  library(boot)
  library(data.table)
  library(doMC)
  
  #Select wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SI" & se$microbiome %in% c("wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
  meth <- t(meth)
  
  #Get actual age of samples
  age <- colData(se)[,"age"]
  
  #Bind age to meth values
  meth_age <- cbind(age,meth)
  
  #Define function that will be run in bootstrapping
  lasso_boot <- function(){
    
    rsamples <- as.data.table(colData(se))[, .SD[sample(.N, replace = TRUE)], by = age][, sample_id]
    
    rmeth_age <- meth[rsamples,]
    
    fit <- cv.glmnet(x = rmeth_age[,-1], y = rmeth_age[,1], alpha = .8)
    
    return(which(coef(fit) > 0))
  }
  
  #Define how many boots and cores will be used for bootstrapping
  boots <- 1000
  cores <- 16
  
  boot_results <- parallel::mclapply(
    1:boots, 
    function(i) suppressMessages(lasso_boot()),
    mc.cores = cores
  )
  
  job::export(c(boot_results))
},title = "Bootstrap select SI aging loci")

#### Bootstrap select SPL aging loci with strain added####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  library(glmnet)
  library(boot)
  library(data.table)
  library(doMC)
  
  #Select wt microbiome samples and perform cov/var filter
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
  se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw") %>% as.matrix()
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get just the wt samples
  se <- se[,se$microbiome %in% c("wt")]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw") %>% as.matrix()
  meth <- t(meth)
  
  #Get actual age of samples
  age <- colData(se)[,"age"]
  
  #Get strain of samples
  strain <- colData(se)[,"strain"] %>% as.factor %>% as.numeric()
  
  #Bind age and strain to meth values
  meth <- cbind(age,strain,meth)
  
  #Define function that will be run in bootstrapping
  lasso_boot <- function(){
    
    rsamples <- as.data.table(colData(se))[, .SD[sample(.N, replace = TRUE)], by = age][, sample_id]
    
    rmeth <- meth[rsamples,]
    
    fit <- cv.glmnet(x = rmeth[,-1], y = rmeth[,1], alpha = .8)
    
    return(which(coef(fit) > 0))
  }
  
  #Define how many boots and cores will be used for bootstrapping
  boots <- 1000
  cores <- 16
  
  boot_results <- parallel::mclapply(
    1:boots, 
    function(i) suppressMessages(lasso_boot()),
    mc.cores = cores
  )
  
  #Convert sites into granges
  sites <- unique(unlist(boot_results))
  granges <- se[sites, ] %>% granges()
  
  #Save results
  dir.create(here("results/rds/15_bootstrap_lasso/"), showWarnings = FALSE)
  saveRDS(boot_results, here("results/rds/15_bootstrap_lasso/boot_results_spl_strain.rds"))
  saveRDS(granges, here("results/rds/15_bootstrap_lasso/granges_spl_strain.rds"))
  job::export(c(boot_results))
},title = "Bootstrap select SPL aging loci add strain to model")


#### Bootstrap select spleen aging loci ####
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  library(glmnet)
  library(boot)
  library(data.table)
  library(doMC)
  
  #Load Summarized Experiment
  se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  
  #obtain the samples Himani used for aging sites
  samples <- read.csv(here("data/metadata/annotated_perm_intestine44_age.tsv"), sep = "\t") %>% 
    colnames() %>% 
    grep(pattern = ".cov.gz", value = TRUE) %>% 
    tools::file_path_sans_ext() %>% 
    tools::file_path_sans_ext()
  
  #filter se for Himani's samples
  se <- se[,samples]
  
  #perform cov/var filter
  se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "all", progressbar = FALSE)
  
  #Filter SE for NA's in meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
  se <- se[!matrixStats::rowAnyNAs(meth),]
  
  #Get resulting meth matrix
  meth <- getMeth(se, type = "raw")
  meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
  meth <- t(meth)
  
  #Get actual age of samples
  age <- colData(se)[,"age"]
  
  #Bind age to meth values
  meth_age <- cbind(age,meth)
  
  #Define function that will be run in bootstrapping
  lasso_boot <- function(){
    
    rsamples <- as.data.table(colData(se))[, .SD[sample(.N, replace = TRUE)], by = age][, sample_id]
    
    rmeth_age <- meth[rsamples,]
    
    fit <- cv.glmnet(x = rmeth_age[,-1], y = rmeth_age[,1], alpha = 1)
    
    return(which(coef(fit) > 0))
  }
  
  #Define how many boots and cores will be used for bootstrapping
  boots <- 1000
  cores <- 8
  
  boot_results <- parallel::mclapply(
    1:boots, 
    function(i) suppressMessages(lasso_boot()),
    mc.cores = cores
  )
  
  #Convert sites into granges
  sites <- unique(unlist(boot_results))
  granges <- se[sites, ] %>% granges()
  
  #Save results
  dir.create(here("results/rds/15_bootstrap_lasso/"), showWarnings = FALSE)
  saveRDS(boot_results, here("results/rds/15_bootstrap_lasso/boot_results_himani.rds"))
  saveRDS(granges, here("results/rds/15_bootstrap_lasso/granges_himani.rds"))
  job::export(c(boot_results))
},title = "Bootstrap select Himani samples aging loci")
