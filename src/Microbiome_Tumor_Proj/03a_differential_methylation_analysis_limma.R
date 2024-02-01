library(here)
source(here::here("scripts", "01_data_preprocessing.R"))
library(iSEEde)

### COMPARE SIDE ###
#Load summarized experiment
pj_normal <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs"), prefix = "pj_normal")

md <- data.frame(colData(pj_normal[,get_Duplicates(pj_normal$mouse_id)]))
md$suborgan_microbiome <- paste(md$suborgan, md$microbiome, sep = "_")
md$mouse_id <- as.factor(md$mouse_id)

#Get methylation values and convert methylation values to mvals
meth <- assay(pj_normal, "Meth") %>% realize_Parallel(workers = 12, nblocks = 6)
mvals <- beta2m(meth[,md$sample_id])

#This design can be used to test paired sample effects
## design <- model.matrix(~0 + suborgan_microbiome + mouse_id, data = md)
design <- model.matrix(~0 + suborgan_microbiome, data = md)
colnames(design) <- gsub(pattern = "suborgan_microbiome", "", colnames(design))

#estimate the correlation between the repeat samples for the same individual
#dupcor <- limma::duplicateCorrelation(mvals, design, block=md$mouse_id)

#Fit mvals to design
fit <- limma::lmFit(mvals, design)

#Create contrast matrix
contrast.matrix <- limma::makeContrasts(
  by_suborgan = (PCOL_consortium + PCOL_spf + PCOL_gf + PCOL_eckp)/4 - (DCOL_consortium + DCOL_spf + DCOL_gf + DCOL_eckp)/4,
  spf_by_suborgan = PCOL_spf - DCOL_spf,
  gf_by_suborgan = PCOL_gf - DCOL_gf,
  eckp_by_suborgan = PCOL_eckp - DCOL_eckp,
  consortium_by_suborgan = PCOL_consortium - DCOL_consortium,
  cons_gf = (PCOL_consortium + DCOL_consortium)/2 - (PCOL_gf + DCOL_gf)/2,
  cons_gf_pcol = PCOL_consortium - PCOL_gf,
  cons_gf_dcol = DCOL_consortium - DCOL_gf,
  spf_gf = (PCOL_spf + DCOL_spf)/2 - (PCOL_gf + DCOL_gf)/2,
  spf_gf_pcol = PCOL_spf - PCOL_gf,
  spf_gf_dcol = DCOL_spf - DCOL_gf,
  eckp_gf = (PCOL_eckp + DCOL_eckp)/2 - (PCOL_gf + DCOL_gf)/2,
  eckp_gf_pcol = PCOL_eckp - PCOL_gf,
  eckp_gf_dcol = DCOL_eckp - DCOL_gf,
  levels=design
)

fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

dm_results_side <- mcsapply(
  colnames(contrast.matrix), 
  function(x){
    limma::topTable(fit2, number = Inf, coef = x, sort.by = "none")
  },
  simplify = FALSE,
  USE.NAMES = TRUE,
  mc.cores = length(colnames(contrast.matrix))
)

dm_results_side <- mcsapply(
  names(dm_results_side),
  FUN = function(x) {
    df <- dm_results_side[[x]]
    df <- dplyr::mutate(df, chr_base = rownames(df))
    CD <- contrastDiff(
      design = design, contrast.matrix = contrast.matrix,
      contrast_name = x, meth = meth[,md$sample_id]
    )
    mcols <- as.data.frame(mcols(pj_normal))
    df <- left_join(df, mcols, by = "chr_base")
    df <- left_join(df, CD, by = "chr_base")
    df$sig <- df$adj.P.Val < 0.1 & abs(df$diff) > .1
    rownames(df) <- df$chr_base
    return(df)
  },
  simplify = FALSE,
  USE.NAMES = TRUE,
  mc.cores = length(dm_results_side)
)

contrasts_side <- data.frame(
  chr_base = mcols(pj_normal)[["chr_base"]],
  row.names = mcols(pj_normal)[["chr_base"]]
)

for(name in names(dm_results_side)){
  df <- dm_results_side[[name]]
  test <- df[df$sig,"chr_base"]
  contrasts_side[[name]] <- contrasts_side$chr_base %in% test
}

### FILTER SIDE AND COMPARE MICROBIOME ###

md <- data.frame(colData(pj_normal[,get_Duplicates(pj_normal$mouse_id)]))
md$suborgan_microbiome <- paste(md$suborgan, md$microbiome, sep = "_")
md$mouse_id <- as.factor(md$mouse_id)

#Get methylation values and convert methylation values to mvals. And filter out sidedness
meth <- assay(pj_normal, "Meth") %>% realize_Parallel(workers = 12, nblocks = 6)

sites <- contrasts_side %>%
  filter(!(gf_by_suborgan | consortium_by_suborgan)) %>%
  pull(chr_base)

meth <- meth[sites,]
mvals <- beta2m(meth[,md$sample_id])

#This design can be used to test paired sample effects
## design <- model.matrix(~0 + suborgan_microbiome + mouse_id, data = md)
design <- model.matrix(~0 + suborgan_microbiome, data = md)
colnames(design) <- gsub(pattern = "suborgan_microbiome", "", colnames(design))

#Fit mvals to design
fit <- limma::lmFit(mvals, design)

#Create contrast matrix
contrast.matrix <- limma::makeContrasts(
  cons_gf_ns = (PCOL_consortium + DCOL_consortium)/2 - (PCOL_gf + DCOL_gf)/2,
  cons_gf_pcol_ns = PCOL_consortium - PCOL_gf,
  cons_gf_dcol_ns = DCOL_consortium - DCOL_gf,
  spf_gf_ns = (PCOL_spf + DCOL_spf)/2 - (PCOL_gf + DCOL_gf)/2,
  spf_gf_pcol_ns = PCOL_spf - PCOL_gf,
  spf_gf_dcol_ns = DCOL_spf - DCOL_gf,
  eckp_gf_ns = (PCOL_eckp + DCOL_eckp)/2 - (PCOL_gf + DCOL_gf)/2,
  eckp_gf_pcol_ns = PCOL_eckp - PCOL_gf,
  eckp_gf_dcol_ns = DCOL_eckp - DCOL_gf,
  levels=design
)

fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

dm_results_mb <- mcsapply(
  colnames(contrast.matrix), 
  function(x){
    limma::topTable(fit2, number = Inf, coef = x, sort.by = "none")
  },
  simplify = FALSE,
  USE.NAMES = TRUE,
  mc.cores = length(colnames(contrast.matrix))
)

dm_results_mb <- mcsapply(
  names(dm_results_mb),
  FUN = function(x) {
    df <- dm_results_mb[[x]]
    df <- dplyr::mutate(df, chr_base = rownames(df))
    CD <- contrastDiff(
      design = design, contrast.matrix = contrast.matrix,
      contrast_name = x, meth = meth[,md$sample_id]
    )
    mcols <- as.data.frame(mcols(pj_normal))
    df <- left_join(df, mcols, by = "chr_base")
    df <- left_join(df, CD, by = "chr_base")
    df$sig <- df$adj.P.Val < 0.1 & abs(df$diff) > .1
    rownames(df) <- df$chr_base
    return(df)
  },
  simplify = FALSE,
  USE.NAMES = TRUE,
  mc.cores = length(dm_results_mb)
)

contrasts_mb <- data.frame(
  chr_base = mcols(pj_normal)[["chr_base"]],
  row.names = mcols(pj_normal)[["chr_base"]]
)

for(name in names(dm_results_mb)){
  df <- dm_results_mb[[name]]
  test <- df[df$sig,"chr_base"]
  contrasts_mb[[name]] <- contrasts_mb$chr_base %in% test
}

#create contrasts dataframe laying out what is significant in each contrast
contrasts <- left_join(contrasts_side, contrasts_mb, by = "chr_base")
mcols(pj_normal) <- DataFrame(left_join(data.frame(mcols(pj_normal)), contrasts, by = "chr_base"))

#embed contrast results in se
dm_results <- c(dm_results_side, dm_results_mb)
isee <- pj_normal
for(i in seq_along(dm_results)){
  isee <- embedContrastResults(dm_results[[i]], isee, names(dm_results)[i], class = "limma")
}

message("saving comparison results...")
HDF5Array::saveHDF5SummarizedExperiment(isee, dir = here::here("data/se/rrbs/"), prefix = "isee_filter_side", replace = FALSE)
saveRDS(contrasts, file = here::here("results/Microbiome_Tumor_Proj/03/contrasts_filter_side.rds"))
saveRDS(dm_results, file = here::here("results/Microbiome_Tumor_Proj/03/dm_results_filter_side.rds"))