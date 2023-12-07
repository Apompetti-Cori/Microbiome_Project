library(here)
source(here::here("scripts", "01_data_preprocessing.R"))
library(plyr)

#Load summarized experiment
pj_normal <- HDF5Array::loadHDF5SummarizedExperiment(here::here("data/se/rrbs"), prefix = "pj_normal")

dm_results <- readRDS(here::here("results/Microbiome_Tumor_Proj/03/dm_results.rds"))

# Make constrast matrix
df <- data.frame(vec = names(dm_results), row.names = names(dm_results))

design <- model.matrix(~0+vec, data = df)
colnames(design) <- gsub(pattern = "vec", "", colnames(design))

contrast.matrix <- limma::makeContrasts(
  contrast1 = consortium_by_suborgan - spf_by_suborgan,
  contrast2 = consortium_by_suborgan - eckp_by_suborgan,
  contrast3 = gf_by_suborgan - consortium_by_suborgan,
  contrast4 = gf_by_suborgan - eckp_by_suborgan,
  contrast5 = gf_by_suborgan - spf_by_suborgan,
  contrast6 = eckp_by_suborgan - spf_by_suborgan,
  contrast7 = gf_by_suborgan - (eckp_by_suborgan + spf_by_suborgan + consortium_by_suborgan),
  contrast8 = by_suborgan - gf_by_suborgan,
  contrast9 = by_suborgan - (gf_by_suborgan + consortium_by_suborgan + spf_by_suborgan + eckp_by_suborgan),
  # contrast20 = eckp_gf - spf_gf, 
  # contrast21 = eckp_gf - consortium_gf,
  # contrast22 = spf_gf - consortium_gf,
  # contrast23 = eckp_gf_pcol - eckp_gf_dcol, 
  # contrast24 = spf_gf_pcol - spf_gf_dcol,
  # contrast25 = consortium_gf_pcol - consortium_gf_dcol,
  # contrast26 = (eckp_gf_pcol + eckp_gf_dcol) - (gf_by_suborgan + consortium_by_suborgan + spf_by_suborgan + eckp_by_suborgan),
  # contrast27 = eckp_gf_pcol - spf_gf_pcol,
  # contrast28 = eckp_gf_pcol - consortium_gf_pcol,
  # contrast29 = eckp_gf_dcol - spf_gf_dcol,
  # contrast30 = eckp_gf_dcol - consortium_gf_dcol,
  # contrast31 = spf_gf_pcol - consortium_gf_pcol,
  # contrast32 = spf_gf_dcol - consortium_gf_dcol,
  levels=design
)

overlap_results <- sapply(
  colnames(contrast.matrix),
  FUN = function(x){
    samples <- contrastSamples(design, contrast.matrix, x)
    
    #full join all dataframes in group1 of the contrast
    df1 <- purrr::reduce(
      dm_results[samples[startsWith(names(samples), "A")]],
      dplyr::inner_join, by = 'chr_base'
    )
    
    #calculate Avgdiff of group1 of the contrast
    df1$Avgdiff <- rowMeans(
      df1[,colnames(df1)[startsWith(colnames(df1), "diff")],drop=FALSE]
    )
    
    #calculate Avgsig of group1 of the contrast. (Only true if one element in group1 are true)
    df1$Avgsig <- rowMeans(
      df1[,colnames(df1)[startsWith(colnames(df1), "sig")],drop=FALSE]
    ) >= (1/length(samples[startsWith(names(samples), "A")]))
    
    #add '.A' suffix to all columns without chr_base
    colnames(df1)[colnames(df1) != "chr_base"] <- paste(colnames(df1)[colnames(df1) != "chr_base"],"A",sep=".")
    
    #full join all dataframes in group2 of the contrast
    df2 <- purrr::reduce(
      dm_results[samples[startsWith(names(samples), "B")]], 
      dplyr::inner_join, by = 'chr_base'
    )
    
    #calculate Avgdiff of group2 of the contrast
    df2$Avgdiff <- rowMeans(
      df2[,colnames(df2)[startsWith(colnames(df2), "diff")],drop=FALSE]
    )
    
    #calculate Avgsig of group2 of the contrast. (Only true if all elements in group1 are true)
    df2$Avgsig <- rowMeans(
      df2[,colnames(df2)[startsWith(colnames(df2), "sig")],drop=FALSE]
    ) >= (1/length(samples[startsWith(names(samples), "B")]))
    
    #add '.B' suffix to all columns without chr_base
    colnames(df2)[colnames(df2) != "chr_base"] <- paste(colnames(df2)[colnames(df2) != "chr_base"],"B",sep=".")
    
    #full join group1 and group2
    df3 <- inner_join(df1, df2, by = "chr_base") %>% 
      dplyr::select(chr_base, Avgsig.A, Avgdiff.A, Avgsig.B, Avgdiff.B) %>%
      mutate(
        color = case_when(
          Avgsig.A & !Avgsig.B ~ "forestgreen",
          !Avgsig.A & Avgsig.B ~ "gold2",
          (Avgsig.A & Avgsig.B) & ((sign(Avgdiff.A) != sign(Avgdiff.B))) ~ "purple",
          (Avgsig.A & Avgsig.B) & ((sign(Avgdiff.A) > 0) & (sign(Avgdiff.B) > 0)) ~ "red",
          (Avgsig.A & Avgsig.B) & ((sign(Avgdiff.A) < 0) & (sign(Avgdiff.B) < 0)) ~ "blue"
        )
      )
    
    
    df4 <- df3 %>% filter(Avgsig.A | Avgsig.B)
    df3[df3$chr_base %noin% df4$chr_base, "color"] <- "gray"
    
    df3 <- df3 %>% mutate(
      facet_var = case_when(
        color == "forestgreen" ~ "A",
        color == "gold2" ~ "B",
        color %in% c("red", "blue", "purple") ~ "AB",
        color == "gray" ~ "Neither"
      )
    )
    
    pearson <- round(cor(df4$Avgdiff.A, df4$Avgdiff.B, method = "pearson"), 2)
    size <- nrow(df4) %>% prettyNum(big.mark = ",")
    maxdiff <- round(max(abs(df3$Avgdiff.A), abs(df3$Avgdiff.B)), 1)
    fisher <- table(df3$Avgsig.A, df3$Avgsig.B) %>% 
      fisher.test(alternative = "two.sided")
    fisherp <- round(fisher$p.value, 2)
    
    p <- df3 %>%
      ggplot(
        aes(x = Avgdiff.A, y = Avgdiff.B, color = color)
      ) +
      geom_point(data = df3[df3$color == "gray", ], alpha = .2, mapping = aes(color = color)) +
      geom_point(data = df3[df3$color %in% c("forestgreen", "gold2"), ], alpha = .4, mapping = aes(color = color)) +
      geom_point(data = df3[df3$color %in% c("red", "blue", "purple"), ], alpha = .4, mapping = aes(color = color)) +
      scale_color_identity() + 
      scale_x_continuous(breaks = seq(-1,1,.2), limits = c(-0.8, 0.8)) + 
      scale_y_continuous(breaks = seq(-1,1,.2), limits = c(-0.8, 0.8)) + 
      labs(
        title = bquote(bold("Overlap of differential methylation")), 
        subtitle = bquote(bold("Overlap size:") ~ .(size) ~ bold("Pearson's coef:") ~ .(pearson) ~ bold("Fisher's p-value:") ~ .(fisherp))
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(face = "bold"),
            axis.text.y = element_text(face = "bold"))
    
    df3$facet_var <- factor(df3$facet_var, levels = c("Neither", "B", "A", "AB"))
    df_counts <- ddply(df3, .(facet_var), summarise, count = length(chr_base))
    
    fp <- ggplot(df3, aes(Avgdiff.A, Avgdiff.B)) + 
      geom_point(alpha = .4, mapping = aes(color=color)) +
      scale_color_identity() +
      theme_classic() +
      facet_wrap(~ facet_var) + 
      scale_x_continuous(breaks = seq(-1,1,.2), limits = c(-0.8, 0.8)) + 
      scale_y_continuous(breaks = seq(-1,1,.2), limits = c(-0.8, 0.8)) + 
      geom_text(data = df_counts, 
                aes(x = Inf, y = Inf, label = paste("n =", count)),
                vjust = "top", hjust = "right", inherit.aes = FALSE, fontface="bold") +
      theme(strip.text = element_text(face = "bold", size = 12),
            axis.text.x = element_text(face = "bold"),
            axis.text.y = element_text(face = "bold"))
    
    mcols <- as.data.frame(mcols(pj_normal))
    df3 <- left_join(df3, mcols, by = "chr_base")
    
    p2 <- df3 %>%
      filter(facet_var == "AB") %>%
      ggplot(mapping = aes(x = color, fill = factor(is_island, levels = c(TRUE, FALSE)))) +
      geom_bar() +
      scale_x_discrete(
        labels = c(
          red = "Red",
          blue = "Blue"
        ),
        limits = c("red", "blue")
      ) +
      scale_fill_manual(name = "CpG Type", values = c("gold", "navy"), labels = c("CpGI", "nonCpGI")) +
      xlab("Color") +
      ylab("Frequency") +
      theme_classic() +
      theme(
        text = element_text(face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 8),
        legend.position = "right"
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 2)),
        shape = guide_legend(override.aes = list(size = 2)),
        fill = guide_legend(override.aes = list(size = 2))
      )
    
    return(list("df" = df3, "plot" = p, "plot_facet" = fp, "ab_plot" = p2))
  },
  simplify = FALSE,
  USE.NAMES = TRUE
)

saveRDS(overlap_results, here::here("results/Microbiome_Tumor_Proj/05/overlap_results.rds"))