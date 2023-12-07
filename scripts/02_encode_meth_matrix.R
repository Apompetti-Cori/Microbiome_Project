job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  
  dt <- data.table::rbindlist(
    parallel::mclapply(bedfiles, data.table::fread, mc.cores = 8), 
    idcol = "sample_id"
  )
  
  dt[, ID := str_c(V1, V2, V2, sep = ".")]
  dt.wider <- data.table::dcast(dt, ID ~ sample_id, value.var = "V11", fill = NA)
  cov_matrix <- data.table::dcast(dt, ID ~ sample_id, value.var = "V10", fill = NA)
  cov_matrix[cov_matrix < 10] <- 0
  #Keep only loci where there is non-zero coverage in at least 2 samples
  keepLoci <- which(rowSums(cov_matrix >= 1) >= 2)
  ref_matrix <- as.matrix(dt.wider, rownames = "ID")
  
  # Remove NA's and convert percents to decimals
  ref_matrix <- (ref_matrix/100)
  ref_matrix.f <- ref_matrix[keepLoci,]
  
  # Remove X chromosome
  ref_matrix.f.x <- ref_matrix.f[!grepl("chrX",rownames(ref_matrix.f)),]
  
  message("saving output")
  saveRDS(ref_matrix.f.x, here("results","rds","ref_matrix.rds"))
  job::export(c(ref_matrix.f.x))
}, title = "Generating ENCODE Meth Matrix")