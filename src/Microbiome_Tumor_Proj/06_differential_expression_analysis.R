source(here::here("scripts", "01_data_preprocessing.R"))
library(edgeR)
library(limma)
library(iSEEde)

rnaseq <- readRDS(here::here("data/se/ProxDistRNAseq/RangedSummarizedExperiment.rds"))
tmp_rnaseq <- readRDS(here::here("data/se/ProxDistRNAseq/SummarizedExperiment.rds"))
rowData(rnaseq) <- rowData(tmp_rnaseq)
coldata <- read.csv(here::here("data/se/ProxDistRNAseq/sample_table_input.csv")) %>% select(sample)
rownames(coldata) <- coldata$sample
coldata$suborgan <- c("PCOL", "PCOL", "DCOL", "PCOL", "DCOL", "DCOL")
coldata$group <- coldata$suborgan
colData(rnaseq) <- DataFrame(coldata)
rm(tmp_rnaseq)

dge <- edgeR::SE2DGEList(subset(rnaseq, isCoding))

dge$gene$gene_length <- NA
dge$genes$gene_length <- max(width(granges(subset(rnaseq, isCoding))))

design <- model.matrix(~ 0 + suborgan, data = dge$samples)
colnames(design) <- gsub(pattern = "suborgan", "", colnames(design))

contrast.matrix <- limma::makeContrasts(
  by_suborgan = PCOL - DCOL,
  levels=design
)

keep <- edgeR::filterByExpr(dge, design = design)

dge <- dge[keep, ,keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

cpms <- edgeR::cpm(dge)

logcpms <- edgeR::cpm(dge, log = TRUE)

dge <- estimateDisp(dge, design = design, robust = TRUE)

message("calculating group cpm/rpkm")
grouplogcpm <- data.frame(edgeR::cpmByGroup(dge, log = TRUE))
colnames(grouplogcpm) <- c("avg_logCPM.DCOL", "avg_logCPM.PCOL")
grouplogcpm$FeatureID <- rownames(grouplogcpm)

grouplogrpkm <- data.frame(edgeR::rpkmByGroup(dge, log = FALSE, gene.length = dge$genes$gene_length))
colnames(grouplogrpkm) <- c("avg_RPKM.DCOL", "avg_RPKM.PCOL")
grouplogrpkm$FeatureID <- rownames(grouplogrpkm)

dge$genes <- left_join(dge$genes, grouplogcpm, by = "FeatureID")
dge$genes <- left_join(dge$genes, grouplogrpkm, by = "FeatureID")

fit <- glmQLFit(dge, design, robust = TRUE)

qlf <- glmQLFTest(fit, contrast = contrast.matrix[,"by_suborgan"])

qlftable <- topTags(qlf, n = Inf, sort.by = "none")
qlftable$table$sig <- (abs(qlftable$table$logFC) > 1 & qlftable$table$FDR < .1)

#rowData(rnaseq)[,"logFC.PvD"] <- NA
#rowData(rnaseq)[rownames(qlftable$table),"logFC.PvD"] <- qlftable$table$logFC

#rowData(rnaseq)[,"FDR.PvD"] <- NA
#rowData(rnaseq)[rownames(qlftable$table),"FDR.PvD"] <- qlftable$table$FDR

rowData(rnaseq)[,"sig.PvD"] <- NA
rowData(rnaseq)[rownames(qlftable),"sig.PvD"] <- qlftable$table$sig

edgeR_res <- iSEEedgeRResults(qlftable)

isee <- embedContrastResults(edgeR_res, rnaseq, "PvD")
metadata(isee)$catch <- NULL

message("saving output")
saveRDS(isee, here::here("data/se/ProxDistRNAseq/isee.rds"))
saveRDS(dge, here::here("results/Microbiome_Tumor_Proj/06/PvD_dge.rds"))