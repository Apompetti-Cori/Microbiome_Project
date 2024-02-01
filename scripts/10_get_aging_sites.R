library(tidyverse) |> suppressWarnings()
library(rtracklayer) |> suppressWarnings()
library(here)

mm9tomm10 <- import.chain(here::here("data", "chainfiles","mm9ToMm10.over.chain"))

#Get aging sites
aging <-read_tsv(here::here('data/metadata/annotated_perm_intestine44_age.tsv')) %>%
  dplyr::rename(empirical_p = empirical.p, spearman=cor) %>%
  mutate(Significance = case_when(((empirical_p < 0.01) & (spearman > 0.5 )) ~ 'Hypermethylated',
                                  ((empirical_p < 0.01) & (spearman < -0.5 )) ~ 'Hypomethylated',
                                  ((empirical_p <= 0.05 & empirical_p  >=0.01)& (spearman > 0.5|spearman < -0.5 )) ~ '2_mod_sig',
                                  TRUE ~'Not significant')) %>% 
  mutate(feature_id= str_replace(chr_base, "_", "-"))

#These are the aging sites as granges
aging_sites <- aging %>% 
  filter(Significance == "Hypomethylated" | Significance == "Hypermethylated") %>% 
  select(chr_base) %>% 
  separate(chr_base, c("seqnames", "start"), "_") %>% 
  mutate("start" = as.integer(start),
         "end" = as.integer(start),
         "width" = 1) %>% 
  makeGRangesFromDataFrame()



seqlevelsStyle(aging_sites) = "UCSC"
  
aging_sites2 <- liftOver(aging_sites, mm9tomm10)
aging_sites2 <- unlist(aging_sites2)
aging_sites$chr_base <- paste0("chr",seq(aging_sites),"_",start(aging_sites))
aging_sites2$chr_base <- paste0("chr",seq(aging_sites2),"_",start(aging_sites2))
genome(aging_sites2) = "mm10"
saveRDS(aging_sites2, file = here::here("results/rds/10/age_sites_mm10.rds"))
saveRDS(aging_sites, file = here::here("results/rds/10/age_sites_mm10.rds"))