source(here::here("scripts", "01_data_preprocessing.R"))
library(glmnet)
library(data.table)
library(boot)
library(doMC)

#Select wt microbiome samples and perform cov/var filter
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
se <- se[, se$organ == "SPL" & se$sample_id != "liv_4402"]
se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "all", progressbar = FALSE)

#Filter SE for NA's in meth matrix
meth <- getMeth(se, type = "raw")
meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
se <- se[!matrixStats::rowAnyNAs(meth),]

#Get resulting meth matrix
meth <- getMeth(se, type = "raw")
meth <- realize_Parallel(meth, cores = 12, nblocks = 6)
meth <- t(meth)

#Get actual age of samples
age <- colData(se)[,"age"]

meth <- cbind(age,meth)

fit2 <- glmnet::cv.glmnet(x = meth[,-1], y = meth[,1], alpha = 1)

se <- se[which(coef(fit) > 0),]

meth <- getMeth(se, type = "raw")
meth <- as.matrix(meth)

pca <- PCAtools::pca(meth, metadata = colData(se))
p <- PCAtools::biplot(pca, colby = "age")
ggplotly(p)
plotly::plot_ly(pca$rotated, x = ~PC1, y = ~PC2, z = ~PC3, 
                color = pca$metadata$age, symbol = pca$metadata$microbiome,
                symbols = c("circle", "square")) 
