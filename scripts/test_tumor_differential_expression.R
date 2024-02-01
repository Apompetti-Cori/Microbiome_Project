se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs/"))
pj_normal <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs/"), prefix = "pj_normal")
rowData(se)$chr_base <- paste(seqnames(se), start(se), sep = "_")
rownames(se) <- rowData(se)$chr_base
se <- se[rownames(pj_normal),]
se <- se[,se$strain == "129svev" & se$treatment == "None" & se$genotype != "WT" & se$microbiome == "spf" & se$mouse_id != "5168"]

md <- data.frame(colData(se[,get_Duplicates(se$mouse_id)]))
md$suborgan_tumor <- paste(md$suborgan, md$sampletype, sep = "_")
md$mouse_id <- as.factor(md$mouse_id)
colnames(design) <- gsub(pattern = "suborgan_tumor", "", colnames(design))

#Get methylation values and convert methylation values to mvals
meth <- assay(se, "Meth") %>% realize_Parallel(workers = 12, nblocks = 6)
mvals <- beta2m(meth[,md$sample_id])

design <- model.matrix(~0 + suborgan_tumor, data = md)

fit <- limma::lmFit(mvals, design)

contrast.matrix <- limma::makeContrasts(
  pcol_by_tumor = PCOL_tumor - PCOL_normal,
  dcol_by_tumor = DCOL_tumor - DCOL_normal,
  by_tumor = (PCOL_tumor + DCOL_tumor)/2 - (PCOL_normal + DCOL_normal)/2,
  levels=design
)

fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

dm_results <- mcsapply(
  colnames(contrast.matrix), 
  function(x){
    df <- limma::topTable(fit2, number = Inf, coef = x, sort.by = "none")
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
  mc.cores = length(colnames(contrast.matrix))
)