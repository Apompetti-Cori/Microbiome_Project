library(here)
source(here::here("scripts", "01_data_preprocessing.R"))

#Load summarized experiment
pj_normal <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs"), prefix = "isee_filter_side")
meth <- assay(pj_normal, "Meth") %>% realize_Parallel(workers = 12, nblocks = 6)

#Load in contrasts
contrasts <- as.data.frame(mcols(pj_normal))

#PCA ALL
pca_meth <- meth
pca_all <- PCAtools::pca(
  pca_meth[!rowAnyNAs(pca_meth),],
  metadata = colData(pj_normal)
)

#PCA NONCGI
sites <- rownames(subset(pj_normal, !is_island))
pca_meth <- meth[sites,]
pca_noncgi <- PCAtools::pca(
  pca_meth[!rowAnyNAs(pca_meth),],
  metadata = colData(pj_normal)
)

#PCA CGI
sites <- rownames(subset(pj_normal, is_island))
pca_meth <- meth[sites,]
pca_cgi <- PCAtools::pca(
  pca_meth[!rowAnyNAs(pca_meth),],
  metadata = colData(pj_normal)
)

#PCA NO SIDE
sites <- contrasts %>% 
  filter(!(gf_by_suborgan | consortium_by_suborgan)) %>%
  pull(chr_base)

pca_meth <- meth[rownames(meth) %in% sites,]
pca_meth <- pca_meth[complete.cases(pca_meth),]

pca_noside <- PCAtools::pca(
  pca_meth,
  metadata = colData(pj_normal)[colnames(pca_meth),]
)

#PCA NO SIDE NONCGI
sites <- contrasts %>% 
  filter(!(gf_by_suborgan | consortium_by_suborgan) & !is_island) %>%
  pull(chr_base)

pca_meth <- meth[rownames(meth) %in% sites,]
pca_meth <- pca_meth[complete.cases(pca_meth),]

pca_noside_noncgi <- PCAtools::pca(
  pca_meth,
  metadata = colData(pj_normal)[colnames(pca_meth),]
)

#PCA NO SIDE CGI
sites <- contrasts %>% 
  filter(!(gf_by_suborgan | consortium_by_suborgan) & is_island) %>%
  pull(chr_base)

pca_meth <- meth[rownames(meth) %in% sites,]

pca_noside_cgi <- PCAtools::pca(
  pca_meth[!rowAnyNAs(pca_meth),],
  metadata = colData(pj_normal)[colnames(pca_meth),]
)

#PCA SIDE
sites <- contrasts %>% 
  filter((gf_by_suborgan | consortium_by_suborgan)) %>%
  pull(chr_base)

pca_meth <- meth[rownames(meth) %in% sites,]

pca_side <- PCAtools::pca(
  pca_meth[!rowAnyNAs(pca_meth),],
  scale = TRUE,
  metadata = colData(pj_normal)[colnames(pca_meth),]
)

#PCA SIDE NONCGI
sites <- contrasts %>% 
  filter((gf_by_suborgan | consortium_by_suborgan) & !is_island) %>%
  pull(chr_base)

pca_meth <- meth[rownames(meth) %in% sites,]

pca_side_noncgi <- PCAtools::pca(
  pca_meth[!rowAnyNAs(pca_meth),],
  scale = TRUE,
  metadata = colData(pj_normal)[colnames(pca_meth),]
)

#PCA SIDE CGI
sites <- contrasts %>% 
  filter((gf_by_suborgan | consortium_by_suborgan) & is_island) %>%
  pull(chr_base)

pca_meth <- meth[rownames(meth) %in% sites,]

pca_side_cgi <- PCAtools::pca(
  pca_meth[!rowAnyNAs(pca_meth),],
  scale = TRUE,
  metadata = colData(pj_normal)[colnames(pca_meth),]
)

pca_results <- list(
  "pca_all" = pca_all,
  "pca_noncgi" = pca_noncgi,
  "pca_cgi" = pca_cgi,
  "pca_noside" = pca_noside,
  "pca_noside_noncgi" = pca_noside_noncgi,
  "pca_noside_cgi" = pca_noside_cgi,
  "pca_side" = pca_side,
  "pca_side_noncgi" = pca_side_noncgi,
  "pca_side_cgi" = pca_side_cgi
)

saveRDS(pca_results, here::here("results/Microbiome_Tumor_Proj/04/pca_results.rds"))