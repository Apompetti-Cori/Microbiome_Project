se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs/"))
pj_normal <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs/"), prefix = "pj_normal")
rowData(se)$chr_base <- paste(seqnames(se), start(se), sep = "_")
rownames(se) <- rowData(se)$chr_base
se <- se[rownames(pj_normal),]
se <- se[,se$strain == "129svev" & se$treatment == "None" & se$genotype != "WT"]

meth <- assay(se, "Meth") %>% realize_Parallel(workers = 12, nblocks = 6)


pca <- PCAtools::pca(meth[!matrixStats::rowAnyNAs(meth),], metadata = colData(se))
