library(here)
library(foreach)
library(iSEEde)
source(here::here("scripts", "01_data_preprocessing.R"))

message("loading bigwigs...")
bigwig_files <- list.files(here::here("data/bigwig"), full.names = TRUE, pattern = "GSM56523[0-9]*\\.hg38\\.bigwig")
beta_files <- list.files(here::here("data/beta"), full.names = TRUE, pattern = "GSM56523[0-9]*\\.hg38\\.beta")
samples <- str_extract(bigwig_files, pattern = "GSM[0-9]*")
md <- data.frame(
  sample_id = samples, 
  tissue = c(rep("PCOL",2), rep("DCOL",3)),
  sample_type = rep("epi",5)
)
rownames(md) <- samples

beta_content_list <- foreach(i=beta_files, j = samples) %do% {
  fname = i
  N <- file.info(fname)$size
  content <- matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)
  colnames(content) <- paste0(j, c(".M", ".Cov"))
  content
}

#load in bigwig files as granges
bigwig_gr_list <- foreach(i=bigwig_files) %do% data.frame(rtracklayer::import.bw(i))
bigwig_gr_list <- foreach(i=bigwig_gr_list, j=samples, beta=beta_content_list) %do% {
  colnames(i)[colnames(i) == "score"] <- paste0(j, ".Meth")
  i[[paste0(j, ".M")]] <- beta[,paste0(j, ".M")]
  i[[paste0(j, ".Cov")]] <- beta[,paste0(j, ".Cov")]
  data.frame2GRanges(i,keepColumns = TRUE)
}

names(bigwig_gr_list) <- samples

merged_granges <- Reduce("merge", bigwig_gr_list)
mcols(merged_granges)$chr_base <- paste(seqnames(merged_granges), start(merged_granges), sep = "_")

rm(bigwig_gr_list, beta_content_list, i, j)

meth <- DelayedArray(mcols(merged_granges[,colnames(mcols(merged_granges))[endsWith(colnames(mcols(merged_granges)), ".Meth")]]))
M <- DelayedArray(mcols(merged_granges[,colnames(mcols(merged_granges))[endsWith(colnames(mcols(merged_granges)), ".M")]]))
Cov <- DelayedArray(mcols(merged_granges[,colnames(mcols(merged_granges))[endsWith(colnames(mcols(merged_granges)), ".Cov")]]))
colnames(meth) <- samples
colnames(M) <- samples
colnames(Cov) <- samples
meth[sign(meth) == -1] <- NA

message("creating unfiltered SE...")
se <- bsseq::BSseq(
  M = M,
  Cov = Cov,
  gr = merged_granges[,c("width", "chr_base")],
  sampleNames = samples
)

colData(se) <- DataFrame(md)
rownames(se) <- mcols(se)$chr_base
assays(se)$Meth <- getMeth(se, type = "raw")

message("saving unfiltered SE...")
HDF5Array::saveHDF5SummarizedExperiment(se, here::here("data/se/GSE186458/"), replace = TRUE)
rm(list = ls())
se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/GSE186458"))

se <- filter_rrbs(
  se, 
  percent_coverage = .66, 
  min_coverage = 20, 
  max_coverage = Inf,
  workers = 40, 
  nblocks = 20
)

message("building annotation GRanges object...")
#Build annotations for hg19 cpgs, genes, & promoters
annotations <- annotatr::build_annotations(
  "hg19", c("hg19_cpgs")
)

#reformat annotations
annotations <- annotations %>%
  data.frame() %>%
  dplyr::select(seqnames, start, end, width, strand, id) %>%
  mutate(id = gsub(pattern = ":[0-9]*", replacement = "", id),
         type = "cpg") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

rownames(mcols(annotations)) <- paste(seqnames(annotations), start(annotations), sep = "_")

#Get gtf file for annotating cpgs
gtf <- rtracklayer::import(here::here("data/metadata/gencode.v44.annotation.gtf.gz"))

transcripts <- subset(gtf, gene_type == "protein_coding" & type == "transcript") %>%
  data.frame() %>% 
  dplyr::select(seqnames, start, end, width, strand, type, gene_id, gene_name) %>%
  dplyr::rename(id = type) %>% 
  dplyr::mutate(type = "transcript") %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

tx_promoters <- promoters(transcripts, upstream = 500, downstream = 500) %>% 
  data.frame %>%
  mutate(id = "promoter",
         type = "transcript") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

genes <- subset(gtf, gene_type == "protein_coding" & type == "gene") %>%
  data.frame() %>% 
  dplyr::select(seqnames, start, end, width, strand, type, gene_id, gene_name) %>%
  dplyr::rename(id = type) %>% 
  dplyr::mutate(type = "gene") %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

#Get promoters of gtf file for annotating cpgs
promoters <- promoters(genes, upstream = 500, downstream = 500) %>% 
  data.frame %>%
  mutate(id = "promoter",
         type = "gene") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

#Merge all granges to be used for annotations
granges_merge <- full_join(data.frame(genes), data.frame(promoters), by = c("seqnames", "start", "end", "width", "strand", "id", "type", "gene_id", "gene_name"))
granges_merge <- full_join(granges_merge, data.frame(transcripts), by = c("seqnames", "start", "end", "width", "strand", "id", "type", "gene_id", "gene_name"))
granges_merge <- full_join(granges_merge, data.frame(tx_promoters), by = c("seqnames", "start", "end", "width", "strand", "id", "type", "gene_id", "gene_name"))
granges_merge <- full_join(granges_merge, data.frame(annotations), by = c("seqnames", "start", "end", "width", "strand", "id", "type"))
granges_merge <- granges_merge %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

message("annotating SummarizedExperiment...")
#Annotate cpgs with merged granges object
annotated <- annotatr::annotate_regions(
  granges(se)[,"chr_base"], 
  granges_merge,
  ignore.strand = TRUE,
  quiet = TRUE
) %>% 
  data.frame() %>%
  dplyr::rename(annot_strand = annot.strand) %>%
  dplyr::select(-annot.width, -annot.start, -annot.end, -annot.seqnames)

#Remove "annot." prefix from columns
colnames(annotated) <- gsub(x = colnames(annotated), pattern = "annot\\.", replacement = "")
annotated$chr.base <- annotated$chr_base

#Make grangeslist from annotated
grl <- makeGRangesListFromDataFrame(annotated, split.field = "chr.base", keep.extra.columns=TRUE)

replacement <- granges(se) %>%
  data.frame() %>%
  filter(chr_base %in% setdiff(rownames(se),names(grl))) %>%
  mutate(chr.base = chr_base,
         id = "no_annot",
         type = "no_annot",
         gene_id = NA,
         gene_name = NA,
         annot_strand = NA) %>%
  select(colnames(annotated))

replacement <- makeGRangesListFromDataFrame(replacement, split.field = "chr.base", keep.extra.columns=TRUE)

grl <- c(grl, replacement)
mcols(grl)$chr_base <- names(grl)
mcols(grl)$is_island <- FALSE
mcols(grl)[subset(unlist(grl), id == "island") %>% names() %>% unique(),]$is_island <- TRUE

grl <- grl[rownames(se)]

message("creating filtered SE...")
se.f <- SummarizedExperiment(
  assays = list("Meth" = assay(se, "Meth"), "M" = assay(se, "M"), "Cov" = assay(se, "Cov")),
  rowData = grl,
  colData = colData(se)
)

message("saving filtered SE...")
HDF5Array::saveHDF5SummarizedExperiment(se.f, here::here("data/se/GSE186458/"), prefix = "filter", replace = FALSE)
rm(list = ls())
se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/GSE186458"), prefix = "filter")

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
dm_results$chr_base <- rownames(dm_results)
CD <- contrastDiff(
  design = design, contrast.matrix = contrast.matrix,
  contrast_name = "by_suborgan", meth = meth
)
dm_results <- left_join(dm_results, CD, by = "chr_base")
rownames(dm_results) <- dm_results$chr_base

message("creating iSEE...")
#embed contrast results in se
isee <- se
isee <- embedContrastResults(dm_results, isee, "by_suborgan", class = "limma")
mcols(isee)$sig <- contrastResults(isee, "by_suborgan")$adj.P.Val < 0.1 & abs(contrastResults(isee, "by_suborgan")$diff) > .1

message("saving iSEE...")
#save results
HDF5Array::saveHDF5SummarizedExperiment(isee, dir = here::here("data/se/GSE186458/"), prefix = "isee", replace = FALSE)