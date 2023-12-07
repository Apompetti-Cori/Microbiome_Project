library(here)
library(foreach)
source(here::here("scripts", "01_data_preprocessing.R"))

bigwig_files <- list.files(here::here("data/bigwig"), full.names = TRUE, pattern = "GSM56523[0-9]*\\.hg38\\.bigwig")
samples <- str_extract(bigwig_files, pattern = "GSM[0-9]*")
md <- data.frame(
  sample_id = samples, 
  tissue = c(rep("PCOL",2), rep("DCOL",3)),
  sample_type = rep("epi",5)
)
rownames(md) <- samples
bigwig_gr_list <- foreach(i=bigwig_files) %do% data.frame(rtracklayer::import.bw(i))
bigwig_gr_list <- foreach(i=bigwig_gr_list, j=samples) %do% {
  colnames(i)[colnames(i) == "score"] <- j
  data.frame2GRanges(i,keepColumns = TRUE)
}
names(bigwig_gr_list) <- samples

merged_granges <- Reduce("merge", bigwig_gr_list)
mcols(merged_granges)$chr_base <- paste(seqnames(merged_granges), start(merged_granges), sep = "_")

rm(bigwig_gr_list, i, j)

meth <- DelayedArray(mcols(merged_granges[,colnames(mcols(merged_granges))[startsWith(colnames(mcols(merged_granges)), "G")]]))
meth[sign(meth) == -1] <- NA

se <- SummarizedExperiment(
  assays = list("Meth" = meth),
  rowData = merged_granges[,c("width", "chr_base")],
  colData = md
)

rownames(se) <- mcols(se)$chr_base

message("saving SE...")
HDF5Array::saveHDF5SummarizedExperiment(se, here::here("data/se/GSE186458/"), replace = TRUE)
rm(list = ls())
se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/GSE186458"))
rmchr <- chrSelectBSseq(se, seqnames = c("chrX", "chrY", "chrM")) %>% granges()
se <- subsetByOverlaps(se, rmchr, invert = TRUE); rm(rmchr)

message("performing DM analysis...")
meth <- assay(se, "Meth") %>% realize_Parallel(workers = 12, nblocks = 6)
meth <- meth[!rowAnyNAs(meth),]
mvals <- beta2m(meth)

design <- model.matrix(~0 + tissue, data = colData(se))
fit <- limma::lmFit(mvals, design)

contrast.matrix <- limma::makeContrasts(
  by_suborgan = (tissuePCOL) - (tissueDCOL),
  levels=design
)

fit2 <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

dm_results <- limma::topTable(fit2, number = Inf, coef = "by_suborgan", sort.by = "none")

#embed contrast results in se
isee <- pj_normal
isee <- embedContrastResults(dm_results, isee, "by_suborgan", class = "limma")

message("saving iSEE...")
#save results
HDF5Array::saveHDF5SummarizedExperiment(isee, dir = here::here("data/se/GSE186458/"), prefix = "isee", replace = FALSE)
