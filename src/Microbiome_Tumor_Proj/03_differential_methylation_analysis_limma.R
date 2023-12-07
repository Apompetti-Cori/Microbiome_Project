library(here)
source(here::here("scripts", "01_data_preprocessing.R"))
library(iSEEde)

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
  levels=design
)

fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

dm_results <- mcsapply(
  colnames(contrast.matrix), 
  function(x){
    limma::topTable(fit2, number = Inf, coef = x, sort.by = "none")
  },
  simplify = FALSE,
  USE.NAMES = TRUE,
  mc.cores = length(colnames(contrast.matrix))
)

dm_results <- mcsapply(
  names(dm_results),
  FUN = function(x) {
    df <- dm_results[[x]]
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
  mc.cores = length(dm_results)
)

contrasts <- data.frame(
  chr_base = mcols(pj_normal)[["chr_base"]],
  row.names = mcols(pj_normal)[["chr_base"]]
)

for(name in names(dm_results)){
  df <- dm_results[[name]]
  test <- df[df$sig,"chr_base"]
  contrasts[[name]] <- contrasts$chr_base %in% test
}

mcols(pj_normal) <- DataFrame(left_join(data.frame(mcols(pj_normal)), contrasts, by = "chr_base"))

contrasts <- data.frame(mcols(pj_normal))

rownames(contrasts) <- contrasts$chr_base

#embed contrast results in se
isee <- pj_normal
for(i in seq_along(dm_results)){
  isee <- embedContrastResults(dm_results[[i]], isee, names(dm_results)[i], class = "limma")
}

HDF5Array::saveHDF5SummarizedExperiment(isee, dir = here::here("data/se/rrbs/"), prefix = "isee", replace = FALSE)
saveRDS(contrasts, file = here::here("results/Microbiome_Tumor_Proj/03/contrasts.rds"))
saveRDS(dm_results, file = here::here("results/Microbiome_Tumor_Proj/03/dm_results.rds"))