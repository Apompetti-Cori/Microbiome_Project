library(here)
source(here::here("scripts/01_data_preprocessing.R"))

#Load SE
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here::here("data/se/rrbs"))

se <- se[,se$genotype %in% c("ApcMin", "ApcMin/IL10KO")]

#Filter project for no treatment
se <- se[,se$treatment == "None"]

message("filtering...")
#Filter project for 20 coverage across 75% of the samples. Ignore tumor samples when filtering.
se <- filter_rrbs(
  se, 
  percent_coverage = .75,
  min_coverage = 20, 
  max_coverage = Inf,
  workers = 40, 
  nblocks = 20,
  ignore_samples = "se$sampletype != 'tumor'"
)

#Filter project for just normal samples.
se <- se[,se$sampletype == "normal"]

#create chr_base column
rowData(se)$chr_base <- ""
rowData(se)$chr_base <- paste(seqnames(se), start(se), sep = "_")
rownames(se) <- rowData(se)$chr_base

message("building annotation GRanges object...")
#Build annotations for mm10 cpgs, genes, & promoters
annotations <- annotatr::build_annotations(
  "mm10", c("mm10_cpgs")
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
gtf <- rtracklayer::import(here::here("data/metadata/gencode.vM25.annotation.gtf.gz"))

transcripts <- subset(gtf, gene_type == "protein_coding" & type == "transcript" & !startsWith(gene_name, "Gm") & !endsWith(gene_name, "Rik")) %>%
  data.frame() %>% 
  dplyr::select(seqnames, start, end, width, strand, type, gene_id, gene_name) %>%
  dplyr::rename(id = type) %>% 
  dplyr::mutate(type = "transcript") %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

tx_promoters <- promoters(transcripts) %>% 
  data.frame %>%
  mutate(id = "promoter",
         type = "transcript") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

genes <- subset(gtf, gene_type == "protein_coding" & type == "gene" & !startsWith(gene_name, "Gm") & !endsWith(gene_name, "Rik")) %>%
  data.frame() %>% 
  dplyr::select(seqnames, start, end, width, strand, type, gene_id, gene_name) %>%
  dplyr::rename(id = type) %>% 
  dplyr::mutate(type = "gene") %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

#Get promoters of gtf file for annotating cpgs
promoters <- promoters(genes) %>% 
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
  granges(se), 
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

#annotated <- annotated %>% dplyr::distinct(chr_base, annot_strand, id, .keep_all = TRUE)

#Make grangeslist from annotated
grl <- makeGRangesListFromDataFrame(annotated, split.field = "chr.base", keep.extra.columns=TRUE)
mcols(grl)$chr_base <- names(grl)
mcols(grl)$is_island <- FALSE
mcols(grl)[subset(unlist(grl), id == "island") %>% names() %>% unique(),]$is_island <- TRUE

ranged_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list("M" = assay(se, "M"), "Cov" = assay(se, "Cov"), "Meth" = assay(se, "Meth")),
  colData = colData(se),
  rowRanges = grl[rownames(se)]
)

message("saving results...")
HDF5Array::saveHDF5SummarizedExperiment(
  ranged_se, 
  dir = here::here("data/se/rrbs"),
  replace = FALSE, 
  prefix = "pj_normal"
)