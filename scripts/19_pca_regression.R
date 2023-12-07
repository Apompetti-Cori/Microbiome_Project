#Select gf/wt microbiome samples and perform cov/var filter
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
se <- se[,se$organ == "SPL" & se$microbiome %in% c("gf","wt") & se$sample_id != "liv_4402"]
granges <- readRDS(here("results/rds/18_lm_across_sites/granges_spl_strain.rds"))
se <- IRanges::subsetByOverlaps(se, granges, invert = FALSE)
meth <- getMeth(se, type = "raw") %>% as.matrix()
pca <- prcomp(t(meth[,se$microbiome == "wt"]))

train <- cbind(pca$x, colData(se[,se$microbiome == "wt"])[,"age",drop=FALSE])
mod <- lm(age~PC1,data = train)

new <- t(meth[,se$microbiome == "gf"])

pred.prcomp <- predict(pca, newdata = new) %>% as.data.frame()

pred.lm <- predict(mod, newdata = pred.prcomp)

pred.lm

plot(pred.lm, colData(se[,se$microbiome == "gf"])[,"age"])
