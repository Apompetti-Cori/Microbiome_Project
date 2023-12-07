source(here::here("scripts", "01_data_preprocessing.R"))
#Load SE with all chromosomes
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
se <- se[,se$sample_id != "liv_4402"]
se <- filter_rrbs(se)
# se <- renameSeqlevels(se, c(MT="M"))
# se <- renameSeqlevels(se, paste("chr",seqlevels(se),sep=""))
# se <- chrSelectBSseq(se, seqnames = c("chrJ02459.1_lambda"))

conv_stats <- read_tsv(here::here("data/metadata/compiled_bs_conversion.tsv"))
conv_stats <- conv_stats %>% rename("sample_id" = "Sample_Name") %>% select(-`library(tidyverse)`)
md <- colData(se) %>% as.data.frame()

x <- data.frame(search = character(), hit = character())
for(pattern in md$sample_id){
  hit <- grep(pattern, conv_stats$sample_id, value = T)
  x <- x %>% add_row(search = pattern, hit = hit)
}

x <- x %>% 
  group_by(search) %>% 
  summarise(sample_id = first(hit))

conv_stats <- merge(x, conv_stats, by = "sample_id") %>% select(-sample_id) %>% rename(sample_id = search)

md <- merge(md, conv_stats, by = "sample_id")
rownames(md) <- md$sample_id

meth <- getMeth(se, type = "raw") %>% as.matrix()
md <- md[colnames(meth),]
meth <- meth[!matrixStats::rowAnyNAs(meth),]
# meth <- meth[matrixStats::rowSums2(is.na(meth))/ncol(meth) < .8,]

pca <- PCAtools::pca(meth, metadata = md)
ggplotly(biplot(pca, colby = "Methylated_nonCpG", shape = "batch"))
ggplotly(pcacorplot(pca, components = 1:10, metavars = c("Methylated_nonCpG", "batch")))
