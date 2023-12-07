


annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

map_ID <- function(metadata, files){
  mapped_meta <- transform(
    metadata,
    filename = unlist(files)[adist(metadata$sample_id, files, fixed = FALSE) %>% 
                               sweep(1,0,"==") %>% 
                               apply(1,function(x) which(x)[1])]
  )
  return(mapped_meta)
}

assign_lapply <- function(named_list){
  sapply(names(named_list), FUN = function(x){
    assign(x, named_list[[x]], envir = .GlobalEnv)
  })
}

get_Duplicates <- function(x, values = FALSE){
  duplicates <- x %in% x[duplicated(x)]
  if(values){
    result <- x[duplicates]
    return(result)
  }
  else{
    return(duplicates)
  }
}

# An mc-version of the sapply function.
mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

contrastDiff <- function(design, contrast.matrix, contrast_name, meth){
  contrast <- contrast.matrix[,contrast_name]
  
  group1 <- rownames(design[rowSums2(design[,contrast > 0, drop = FALSE], useNames = FALSE) == 1,])
  group2 <- rownames(design[rowSums2(design[,contrast < 0, drop = FALSE], useNames = FALSE) == 1,])
  
  mean1 <- rowMeans(meth[,group1], na.rm = TRUE)
  mean2 <- rowMeans(meth[,group2], na.rm = TRUE)
  
  x <- data.frame("group1" = mean1, "group2" = mean2,
                  "diff" = mean1 - mean2, 
                  "chr_base" = rownames(meth), 
                  row.names = rownames(meth))
  
  x$status[x$diff > 0] <- "hyper"
  x$status[x$diff < 0] <- "hypo"
  
  return(x)
}

contrastSamples<- function(design, contrast.matrix, contrast_name, separate = FALSE){
  contrast <- contrast.matrix[,contrast_name, drop = FALSE]
  
  group1 <- rownames(design[rowSums2(design[,contrast > 0, drop = FALSE], useNames = FALSE) == 1,,drop=FALSE])
  group2 <- rownames(design[rowSums2(design[,contrast < 0, drop = FALSE], useNames = FALSE) == 1,,drop=FALSE])
  
  return(c("A" = group1, "B" = group2))
}

init_here <- function() {
  `%noin%` = Negate(`%in%`)
  files <- dir( all.files = T )
  while ( ".here" %noin% files & getwd()!="/" ) {
    setwd("..")
    files <- dir( all.files = T )
  }
  i_am(".here")
}

beta2m <- function(x, offset = 1) {
  log2( (x + offset) / (1 - x + offset) )
}

"%noin%" <- Negate("%in%")

se_structure <- function(md){
  md <- as.data.frame(md)
  for(o in unique(md$organ)){
    md_organ <- md %>% filter(organ == o)
    for(b in unique(md_organ$batch)){
      md_batch <- md_organ %>% filter(batch == b)
      for(s in unique(md_batch$strain)){
        cat(paste(o, b, s),"\n")
        md_strain <- md_batch %>% filter(strain == s)
        age <- md_strain$age
        microbiome <- md_strain$microbiome
        print(table(age, microbiome))
      }
    }
  }
}

#### Delayed Matrix Helpers ####

##### Block functions #####
block.impute.row_mean <- function(block){
  row_means <- DelayedMatrixStats::rowMeans2(block, na.rm = TRUE)
  for(i in 1:nrow(block)){
    nas <- is.na(block[i,])
    if(any(nas)){
      block[i,nas] <- row_means[i]
    }
  }
  return(block)
}

block.row_naSum <- function(block){
  sumnas = 0
  for(i in 1:nrow(block)){
    nas <- is.na(block[i,])
    if(any(nas)){
      sumnas = sumnas + 1
    }
  }
  return(sumnas)
}

block.missingness <- function(delayedmat, progressbar = FALSE){
  nas <- blockApply(delayedmat, grid = rowAutoGrid(delayedmat, nrow = 10000),
                    function(block){sum(is.na(block))},
                    BPPARAM = BiocParallel::MulticoreParam(workers = 10, progressbar = progressbar))
  nas <- sum(unlist(nas))
  return(nas/length(meth_raw_filter_75))
}

block.rowSum <- function(delayedmat, progressbar = FALSE, grid, na.rm = TRUE){
  rowsumblock <- blockApply(delayedmat, grid = grid,
                            function(block){rowSums(block, na.rm = na.rm)},
                            BPPARAM = BiocParallel::MulticoreParam(workers = 10, progressbar = progressbar))
  return(do.call(c, rowsumblock))
}

missingness.delayed <- function(delayedmat){
  sum(is.na(delayedmat))/length(delayedmat) * 100
}

##### Block Apply Functions #####
complete.cases.sum.delayed <- function(delayedmat, dim, progressbar = FALSE){
  if(dim == "row"){
    rowNaSum <- blockApply(delayedmat, grid = rowAutoGrid(delayedmat, nrow = 10000),
               function(block){block.row_naSum(block)},
               BPPARAM = BiocParallel::MulticoreParam(workers = 10, progressbar = progressbar))
    rowNaSum <- sum(unlist(rowNaSum))
    
    return(nrow(delayedmat) - rowNaSum)
  }
}

rowSums_Parallel <- function(DelayedMatrix, workers = 1, nblocks = 1, na.rm = FALSE, progressbar = FALSE){
  rowSums <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks), function(block){DelayedMatrixStats::rowSums2(block, na.rm = na.rm)}, 
                        BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = progressbar))
  return(unlist(rowSums))
}

rowstdev_Parallel <- function(DelayedMatrix, workers = 1, nblocks = 1, progressbar = FALSE){
  row_stdev <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks), 
                        function(block){apply(block, 1, sd, na.rm = TRUE)}, 
                        BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = progressbar))
  return(unlist(row_stdev))
}

rowMeans_Parallel <- function(DelayedMatrix, workers = 1, nblocks = 1, na.rm = FALSE, progressbar = FALSE){
  rowSums <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks), 
                        function(block){DelayedMatrixStats::rowMeans2(block, na.rm = na.rm)}, 
                        BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = progressbar))
  return(unlist(rowSums))
}

rowVars_Parallel <- function(DelayedMatrix, workers = 1, nblocks = 1, na.rm = FALSE, progressbar = FALSE){
  rowVars <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks), 
                        function(block){DelayedMatrixStats::rowVars(block, na.rm = na.rm)}, 
                        BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = progressbar))
  return(unlist(rowVars))
}

rowSds_Parallel <- function(DelayedMatrix, workers = 1, nblocks = 1, na.rm = FALSE, progressbar = FALSE){
  rowSds <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks), 
                        function(block){DelayedMatrixStats::rowSds(block, na.rm = na.rm)}, 
                        BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = progressbar))
  return(unlist(rowSds))
}

rowAnyNA_Parallel <- function(DelayedMatrix, workers = 1, nblocks = 1, progressbar = FALSE){
  naIndices <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks),
                          function(block){DelayedMatrixStats::rowAnyNAs(block, useNames = TRUE)},
                          BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = progressbar))
  return(unlist(naIndices))
}

realize_Parallel <- function(DelayedMatrix, workers = 1, nblocks = 1, progressbar = FALSE){
  chunks <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks),
                       function(block){as.matrix(block)},
                       BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = progressbar))
  
  return(do.call(rbind, chunks))
}


#### SE Functions ####
annotate_se <- function(se){
  library(annotatr)
  
  rowData(se)$chr_base <- ""
  rowData(se)$chr_base <- paste(seqnames(se), start(se), sep = "_")
  rownames(se) <- rowData(se)$chr_base
  
  annotations <- annotatr::build_annotations(
    "mm10", c("mm10_genes_promoters", "mm10_cpgs")
  )
  
  annotated <- annotate_regions(
    regions = granges(se),
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE
  ) %>% data.frame()
  
  message("Annotating gene ID's/symbols...")
  #annotate pos strand gene_id
  annotated_transformed <- annotated %>%
    group_by(chr_base) %>%
    distinct(chr_base, annot.gene_id, .keep_all = TRUE)
  
  test <- annotated_transformed %>% 
    filter(annot.strand == "+")
  rowData(se)$annot.gene_id_pos <- ""
  rowData(se)[test$chr_base,]$annot.gene_id_pos <- test$annot.gene_id
  
  #annotate neg strand gene_id
  test <- annotated_transformed %>% 
    filter(annot.strand == "-")
  rowData(se)$annot.gene_id_neg <- ""
  rowData(se)[test$chr_base,]$annot.gene_id_neg <- test$annot.gene_id
  
  #annotate pos strand symbol
  annotated_transformed <- annotated %>%
    group_by(chr_base) %>%
    distinct(chr_base, annot.symbol, .keep_all = TRUE)
  
  test <- annotated_transformed %>% 
    filter(annot.strand == "+")
  rowData(se)$annot.symbol_pos <- ""
  rowData(se)[test$chr_base,]$annot.symbol_pos <- test$annot.symbol
  
  #annotate neg strand symbol
  test <- annotated_transformed %>% 
    filter(annot.strand == "-")
  rowData(se)$annot.symbol_neg <- ""
  rowData(se)[test$chr_base,]$annot.symbol_neg <- test$annot.symbol
  
  message("Annotating promoters...")
  #annotate cpgs for promoters
  test <- annotated %>% 
    filter(str_detect(annot.id, "promoter")) %>% 
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_promoter <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  #annotate pos strand for promoters
  test <- annotated %>% 
    filter(str_detect(annot.id, "promoter") & annot.strand == "+") %>% 
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_promoter_pos <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  #annotate neg strand for promoters
  test <- annotated %>% 
    filter(str_detect(annot.id, "promoter") & annot.strand == "-") %>% 
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_promoter_neg <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  message("Annotating islands...")
  #annotate cpgs for cpg-islands
  test <- annotated %>% 
    filter(str_detect(annot.id, "island")) %>% 
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_island <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  message("Annotating shores...")
  #annotate cpgs for cpg-shores
  test <- annotated %>% 
    filter(str_detect(annot.id, "shore")) %>% 
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_shore <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  message("Annotating shelves...")
  #annotate cpgs for cpg-shelves
  test <- annotated %>% 
    filter(str_detect(annot.id, "shelf")) %>% 
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_shelf <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  message("Annotating interCGI's...")
  #annotate cpgs for interCGI
  test <- annotated %>% 
    filter(str_detect(annot.id, "inter")) %>% 
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_inter <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  message("Annotating for rep elements...")
  rep_elems <- rtracklayer::import(here::here("data/bedfiles/mm10.fa.out.bed"), format = "bed")
  LINE <- subset(rep_elems, grepl("LINE", rep_elems$name))
  SINE <- subset(rep_elems, grepl("SINE", rep_elems$name))
  LTR <-  subset(rep_elems, grepl("LTR", rep_elems$name))
  DNA_trans <- subset(rep_elems, grepl("DNA", rep_elems$name))
  
  #annotate for any rep element
  overlap <- findOverlaps(rep_elems, granges(se))
  test <- rownames(se)[subjectHits(overlap)]
  rowData(se)$is_rep_elem <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  #annotate for LINE type rep elements
  overlap <- findOverlaps(LINE, granges(se))
  test <- rownames(se)[subjectHits(overlap)]
  rowData(se)$is_LINE <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  #annotate for SINE type rep elements
  overlap <- findOverlaps(SINE, granges(se))
  test <- rownames(se)[subjectHits(overlap)]
  rowData(se)$is_SINE <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  #annotate for LTR type rep elements
  overlap <- findOverlaps(LTR, granges(se))
  test <- rownames(se)[subjectHits(overlap)]
  rowData(se)$is_LTR <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  #annotate for DNA_trans type rep elements
  overlap <- findOverlaps(DNA_trans, granges(se))
  test <- rownames(se)[subjectHits(overlap)]
  rowData(se)$is_DNA_trans <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  message("Annotating for coding regions...")
  #annotate for protein coding regions
  gtf <- rtracklayer::import(here::here("data/metadata/Mus_musculus.GRCm38.102.gtf"))
  gtf <- subset(gtf, type == "gene" & gene_biotype == "protein_coding")
  seqlevels(gtf) <- paste0("chr",seqlevels(gtf))
  
  annotated <- annotatr::annotate_regions(
    regions = granges(se),
    annotations = gtf,
    ignore.strand = TRUE,
    quiet = FALSE
  ) %>% 
    data.frame()
  
  #annotate pos/neg strand for coding region
  test <- annotated %>%
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_coding <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  #annotate pos strand for coding region
  test <- annotated %>% 
    filter(annot.strand == "+") %>%
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_coding_pos <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )
  
  #annotate neg strand for coding region
  test <- annotated %>% 
    filter(annot.strand == "-") %>%
    pull(chr_base) %>% 
    unique()
  rowData(se)$is_coding_neg <- ifelse(
    rowData(se)$chr_base %in% test, TRUE, FALSE
  )

  return(se)
}

filter_rrbs <- function(se, percent_coverage = .66, min_coverage = 10,
                        max_coverage = 500,
                        workers = 5, nblocks = 100, 
                        progressbar = FALSE, ignore_samples = "all()",
                        cpg_select = "allcg"){
  
  message("keeping standard chromosomes...")
  se <- keepStandardChromosomes(se, pruning.mode="coarse")
  
  message("filtering out X/Y/M chromosomes...")
  rmchr <- chrSelectBSseq(se, seqnames = c("chrX","chrY","chrM")) %>% granges()
  se <- subsetByOverlaps(se, rmchr, invert = TRUE)
  
  #Select "All CpGs", "CpG Islands", or "Non CpG Islands"
  if(cpg_select == "cgi"){
    se <- subset(se, is_island == TRUE)
  }
  if(cpg_select == "noncgi"){
    se <-  subset(se, is_island == FALSE)
  }
  if(cpg_select == "allcg"){
    se <- se
  }
  
  cov <- bsseq::getCoverage(se[,eval(parse(text=ignore_samples))], type = "Cov")
  
  message("filtering loci for coverage...")
  message("filtering on ", ncol(se[,eval(parse(text=ignore_samples))]), " samples...")
  keepLoci <- which(rowSums_Parallel(cov >= min_coverage & cov <= max_coverage, workers = workers, nblocks = nblocks, na.rm = TRUE, progressbar = progressbar) >= round(ncol(cov)*percent_coverage))
  
  se <- se[keepLoci,]
  
  # meth <- bsseq::getMeth(se[,eval(parse(text=ignore_samples))], type = "raw")
  
  # message("filtering loci for variance...")
  # message("filtering on ", ncol(se[,eval(parse(text=ignore_samples))]), " samples...")
  # var <- rowVars_Parallel(meth, workers = workers, nblocks = nblocks, na.rm = TRUE, progressbar = progressbar)
  # 
  # se <- se[var != 0,]
  
  return(se)
}

#### PCA Corplot Function ####

pcacorplot <- function(pca, components, metavars){
  
  pc_scores <- cbind(pca$rotated, pca$metadata[,c(metavars), drop = FALSE])
  
  #convert character columns to factors
  pc_scores[sapply(pc_scores, is.character)] <- lapply(pc_scores[sapply(pc_scores, is.character)], as.factor)
  #convert factor columns to numeric for correlation
  pc_scores[sapply(pc_scores, is.factor)] <- lapply(pc_scores[sapply(pc_scores, is.factor)], as.numeric)
  
  corr<- Hmisc::rcorr(as.matrix(pc_scores), type = "pearson")
  cor.mat <- reshape2::melt(corr$r[components,c(metavars),drop = FALSE]^2, 
                            na.rm = TRUE, as.is = FALSE)
  p.mat <- reshape2::melt(corr$P[components,c(metavars),drop = FALSE])
  
  plot.mat <- merge(cor.mat, p.mat, by = c("Var1", "Var2")) %>%
    rename("r" = "value.x", "p" = "value.y") %>%
    mutate(r = round(r, digits = 2),
           label = case_when(p < 0.001 ~ "***",
                             p < 0.01 ~ "**",
                             p < 0.05 ~ "*",
                             TRUE ~ "")) %>% 
    mutate(label = paste(r, label, sep = ""))
  
  label = plot.mat$label
  
  p <- ggplot2::ggplot(data = plot.mat, 
                       mapping = ggplot2::aes_string(x = "Var1", y = "Var2", fill = "r")) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::scale_fill_gradient2(low = "snow", high = "darkgreen", 
                                  mid = "gold", midpoint = .5, limit = c(0, 1), space = "Lab", name = "Corr") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                       vjust = 1, size = 9, hjust = 1), 
                   axis.text.y = ggplot2::element_text(size = 9)) +
    ggplot2::coord_fixed() +
    ggplot2::geom_text(mapping = ggplot2::aes_string(x = "Var1", 
                                                     y = "Var2"), label = label, color = "black", size = 3.5) +
    ggplot2::guides(size = "none") + 
    labs(x = NULL, y = NULL) +
    theme_minimal()
  
  return(p)
}


#### 3D PCA Plot ####

pca3d <- function(pca, perc_var, shape, color){
  plot_ly(data = pca, x = ~PC1, y = ~PC2, z = ~PC3,
          width = 1000, height = 700,
          color = color) %>%
    layout(scene = list(
      xaxis = list(title = paste("PC1", " ", round(perc_var[1],2), "%", sep = "")),
      yaxis = list(title = paste("PC2", " ", round(perc_var[2],2), "%", sep = "")),
      zaxis = list(title = paste("PC3", " ", round(perc_var[3],2), "%", sep = ""))
    ))  %>%
    add_markers(symbol =  shape, 
                symbols = c( "circle", "square", "diamond", "cross", "circle-open", "square-open", "diamond-open", "x"))
}


#### PCAtools to prcomp helper function ####

prcompPCAtools <- function(pca = NULL){
  prcomp <- list(sdev = pca$sdev,
                 rotation = data.matrix(pca$loadings),
                 x = data.matrix(pca$rotated),
                 center = TRUE, scale = TRUE)
  
  class(prcomp) <- 'prcomp'
  
  return(prcomp)
}