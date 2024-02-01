library(here)
library(foreach)
library(doParallel)
library(eulerr)
registerDoParallel(cores = 8)
source(here::here("scripts", "01_data_preprocessing.R"))

#Load summarized experiment
se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs"), prefix = "isee_filter_side")
contrasts <- data.frame(mcols(se)[,!startsWith(colnames(mcols(se)), "iSEE")])

#Declare what kind of comparisons you wish to make
comparisons <- list(
  eckp_spf_consortium = c("eckp_gf_ns", "spf_gf_ns", "cons_gf_ns"),
  eckp_spf_consortium_col = c("eckp_gf_pcol_ns", "spf_gf_pcol_ns", "cons_gf_pcol_ns", "eckp_gf_dcol_ns", "spf_gf_dcol_ns", "cons_gf_dcol_ns"),
  eckp_spf_consortium_dcol = c("eckp_gf_dcol_ns", "spf_gf_dcol_ns", "cons_gf_dcol_ns"),
  eckp_spf_consortium_pcol = c("eckp_gf_pcol_ns", "spf_gf_pcol_ns", "cons_gf_pcol_ns"),
  eckp_pcol_dcol = c("eckp_gf_pcol_ns", "eckp_gf_dcol_ns"),
  spf_pcol_dcol = c("spf_gf_pcol_ns", "spf_gf_dcol_ns"),
  consortium_pcol_dcol = c("cons_gf_pcol_ns", "cons_gf_dcol_ns")
)


p_list <- foreach(comp = names(comparisons)) %do% {
  message(paste0("running comparison: ", comp,"..."))
  type <- length(comparisons[[comp]])
  message(paste0("\tcomparison is of length: ", type))
  if(type == 2){
    euler_list <- list(
      contrasts[contrasts[,comparisons[[comp]][[1]]],"chr_base"],
      contrasts[contrasts[,comparisons[[comp]][[2]]],"chr_base"]
    )
    
    names(euler_list) <- comparisons[[comp]]
    names(euler_list) <- gsub("_ns", "", names(euler_list))
    names(euler_list) <- gsub("_", " ", names(euler_list))
    names(euler_list) <- toupper(names(euler_list))
    p <- plot(euler(euler_list, shape = "ellipse"), quantities = TRUE)
    message(paste0("\tcomparison ", comp, " done..."))
  }
  
  if(type == 3){
    euler_list <- list(
      contrasts[contrasts[,comparisons[[comp]][[1]]],"chr_base"],
      contrasts[contrasts[,comparisons[[comp]][[2]]],"chr_base"],
      contrasts[contrasts[,comparisons[[comp]][[3]]],"chr_base"]
    )
    
    names(euler_list) <- comparisons[[comp]]
    names(euler_list) <- gsub("_ns", "", names(euler_list))
    names(euler_list) <- gsub("_", " ", names(euler_list))
    names(euler_list) <- toupper(names(euler_list))
    p <- plot(euler(euler_list, shape = "ellipse"), quantities = TRUE)
    message(paste0("\tcomparison ", comp, " done..."))
  }
  
  if(type == 6){
    euler_list <- list(
      contrasts[contrasts[,comparisons[[comp]][[1]]],"chr_base"],
      contrasts[contrasts[,comparisons[[comp]][[2]]],"chr_base"],
      contrasts[contrasts[,comparisons[[comp]][[3]]],"chr_base"],
      contrasts[contrasts[,comparisons[[comp]][[4]]],"chr_base"],
      contrasts[contrasts[,comparisons[[comp]][[5]]],"chr_base"],
      contrasts[contrasts[,comparisons[[comp]][[6]]],"chr_base"]
    )  
    
    names(euler_list) <- comparisons[[comp]]
    names(euler_list) <- gsub("_ns", "", names(euler_list))
    names(euler_list) <- gsub("_", " ", names(euler_list))
    names(euler_list) <- toupper(names(euler_list))
    p <- plot(euler(euler_list, shape = "ellipse"), quantities = TRUE)
    message(paste0("\tcomparison ", comp, " done..."))
  }
  
  return(p)
}
names(p_list) <- names(comparisons)

message("saving output...")
saveRDS(p_list, file = here::here("results/Microbiome_Tumor_Proj/05/euler_plots.rds"))