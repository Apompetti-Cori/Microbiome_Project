source(here::here("scripts", "01_data_preprocessing.R"))
HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))

md <- read.csv(here::here("data/metadata/rrbs_metadata.csv"), row.names = 1)
colData <- md %>% filter(sampletype != "Organoid")
covfiles <- paste(here::here("data/covfiles/"), colData %>% pull(filename), sep = "")

se <- bsseq::read.bismark(
  files = covfiles,
  colData = colData,
  verbose = TRUE,
  BACKEND = "HDF5Array",
  dir = here("data/se/rrbs"),
  replace = TRUE,
  BPPARAM = BiocParallel::MulticoreParam(workers = 60, progressbar = TRUE)
)

message("converting seqnames")
se <- renameSeqlevels(se, c(MT="M"))
se <- renameSeqlevels(se, paste("chr",seqlevels(se),sep=""))

assays(se)$Meth <- getMeth(se, type = "raw")

message("saving changes")
HDF5Array::quickResaveHDF5SummarizedExperiment(
  x = se
)
