library(here)
source(here::here("scripts", "01_data_preprocessing.R"))

aging_samples <-read_tsv(here::here('data/metadata/annotated_perm_intestine44_age.tsv'), show_col_types = FALSE)
aging_samples <- tools::file_path_sans_ext(colnames(aging_samples)[endsWith(colnames(aging_samples), ".cov.gz")], compression = TRUE)

### RUN PERM COR ON HIMANI's SAMPLES WITH MY CPGS ###
se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs"))
pj_normal <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs/"),prefix = "pj_normal")

se <- subsetByOverlaps(se[,se$sample_id %in% aging_samples], pj_normal)
rownames(se) <- rownames(pj_normal)

meth <- assay(se, "Meth") %>% 
  realize_Parallel(
    workers = 12, 
    nblocks = 6
  )

rownames(meth) <- rownames(se)

#keep loci that have at least 2 values
keepLoci <- (rowSums(is.na(meth)) < 42)

meth <- meth[keepLoci,]

age_cor <- scrime::rowCors(meth, colData(se)$age)
age_perm_cor <- coriell::permutation_correlation_test(
  X = meth,
  y = colData(se)$age,
  n_core = 12,
  n_perm = 1000,
  method = "spearman"
)

### RUN PERM COR ON HIMANI's SAMPLES WITH HER CPGS ###
se <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs"))
se <- se[,se$sample_id %in% aging_samples]

se <- filter_rrbs(se, 
                  percent_coverage = .75, 
                  min_coverage = 20)

meth <- assay(se, "Meth") %>% 
  realize_Parallel(
    workers = 12, 
    nblocks = 6
  )

rownames(meth) <- rownames(se)

age_perm_cor_2 <- coriell::permutation_correlation_test(
  X = meth,
  y = colData(se)$age,
  n_core = 12,
  n_perm = 1000,
  method = "spearman"
)

age_perm_cor_2 <- age_perm_cor_2 %>%
  mutate(
    significance = case_when(
      ((empirical.p < 0.01) & (cor > 0.5 )) ~ "Hypermethylated",
      ((empirical.p < 0.01) & (cor < -0.5 )) ~ "Hypomethylated",
      ((empirical.p <= 0.05 & empirical.p  >=0.01)& (cor > 0.5 | cor < -0.5 )) ~ "2_mod_sig",
      TRUE ~"Not significant"
    )
  )

message("saving output...")
saveRDS(age_cor, here::here("results/Microbiome_Tumor_Proj/08/age_cor.rds"))
saveRDS(age_perm_cor, here::here("results/Microbiome_Tumor_Proj/08/age_perm_cor_my_sites.rds"))
saveRDS(age_perm_cor_2, here::here("results/Microbiome_Tumor_Proj/08/age_perm_cor_hv_sites.rds"))