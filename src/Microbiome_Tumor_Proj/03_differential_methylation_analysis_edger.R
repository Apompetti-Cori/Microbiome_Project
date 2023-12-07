library(here)
source(here::here("scripts", "01_data_preprocessing.R"))
library(iSEEde)
library(edgeR)
library(gdata)

#Load summarized experiment
pj_normal <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs"), prefix = "pj_normal")

md <- data.frame(colData(pj_normal[,get_Duplicates(pj_normal$mouse_id)])) %>% select(sample_id, suborgan, microbiome, mouse_id)
md$suborgan_microbiome <- paste(md$suborgan, md$microbiome, sep = "_")
md$mouse_id <- as.factor(md$mouse_id)

md1 <- md
md2 <- md
rownames(md1) <- paste(rownames(md1), "Me", sep = ".")
rownames(md2) <- paste(rownames(md2), "Un", sep = ".")

samples <- interleave(md1, md2)

message("calculating counts and libsizes...")
#Get counts in Methylated context and counts in Unmethylated context
M <- assay(pj_normal[,md$sample_id], "M")  %>% realize_Parallel(workers = 12, nblocks = 6)
Cov <- assay(pj_normal[,md$sample_id], "Cov")  %>% realize_Parallel(workers = 12, nblocks = 6)
U <- Cov - M; rm(Cov)
libsizes <- (colSums(M) + colSums(U))/2
colnames(M) <- paste(colnames(M), "Me", sep = ".")
colnames(U) <- paste(colnames(U), "Un", sep = ".")

#Interleve the columns of the two matrices so it is formatted as "Sample.Me, Sample.Un"
counts <- t(interleave(t(M), t(U)))

message("constructing dge...")
genes <- data.frame(
  seqnames = as.character(unique(seqnames(pj_normal))),
  start = as.numeric(unique(start(pj_normal))),
  end = as.numeric(unique(end(pj_normal))),
  chr_base = as.character(rowData(pj_normal)$chr_base)
)

#construct dgelist
dge <- DGEList(counts, lib.size = rep(libsizes, each = 2), samples = samples, group = samples$suborgan_microbiome, genes = genes)

message("constructing design...")
#Create design matrix
design <- modelMatrixMeth(~ 0 + md$suborgan_microbiome)
colnames(design) <- gsub(pattern = "md\\$suborgan_microbiome", "", colnames(design))


message("constructing contrasts...")
#Create contrast matrix
contrast.matrix <- limma::makeContrasts(
  by_suborgan = (PCOL_consortium + PCOL_spf + PCOL_gf + PCOL_eckp)/4 - (DCOL_consortium + DCOL_spf + DCOL_gf + DCOL_eckp)/4,
  spf_by_suborgan = PCOL_spf - DCOL_spf,
  gf_by_suborgan = PCOL_gf - DCOL_gf,
  eckp_by_suborgan = PCOL_eckp - DCOL_eckp,
  consortium_by_suborgan = PCOL_consortium - DCOL_consortium,
  levels=design
)

message("estimating dispersion factors...")
#estimate dispersion factors
dge <- estimateGLMCommonDisp(dge, design, subset = 10000)

message("fitting to design...")
#Fit counts to design
fit <- glmFit(dge, design)

message("performing contrasts...")
lrt <- glmLRT(fit, contrast=contrast.matrix)

message("saving output...")
saveRDS(dge, here::here("results/Microbiome_Tumor_Proj/03/dge.rds"))
saveRDS(fit, here::here("results/Microbiome_Tumor_Proj/03/fit.rds"))
saveRDS(lrt, here::here("results/Microbiome_Tumor_Proj/03/lrt.rds"))