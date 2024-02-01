library(here)
library(foreach)
source(here::here("scripts", "01_data_preprocessing.R"))
#Load summarized experiment
se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs"), prefix = "isee_filter_side")
contrasts <- data.frame(mcols(se)[,!startsWith(colnames(mcols(se)), "iSEE")])

#Create strings of contrasts you'd wish to run through GO enrichment
contrasts_of_interest <- colnames(contrasts) %>% grep("by|gf", ., value = TRUE)

go_results <- foreach(con = contrasts_of_interest) %do% {
  se_subset <- subset(se, eval(parse(text = con)))
  
  if(grepl("suborgan", con)){
    sites <- data.frame(contrastResults(se_subset, con)) %>% pull(chr_base)
    genes <- data.frame(unlist(granges(se_subset[sites,]))) %>%
      filter(id == "promoter") %>%
      pull("gene_id") %>%
      gsub("\\.[0-9]+","",.)
    
    go_result <- clusterProfiler::enrichGO(
      gene = unique(genes),
      OrgDb = "org.Mm.eg.db",
      ont = "BP",
      keyType = "ENSEMBL",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "fdr"
    )
  }
  
  else{
    sites <- data.frame(contrastResults(se_subset, con)) %>% filter(status == "hyper") %>% pull(chr_base)
    genes <- data.frame(unlist(granges(se_subset[sites,]))) %>%
      filter(id == "promoter") %>%
      pull("gene_id") %>%
      gsub("\\.[0-9]+","",.)
    
    go_result <- clusterProfiler::enrichGO(
      gene = unique(genes),
      OrgDb = "org.Mm.eg.db",
      ont = "BP",
      keyType = "ENSEMBL",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "fdr"
    )
  }
  
  return(go_result)
}
names(go_results) <- contrasts_of_interest

message("saving results...")
saveRDS(go_results, file = here::here("results/Microbiome_Tumor_Proj/09/go_results.rds"))