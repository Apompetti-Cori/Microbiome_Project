source(here::here("scripts", "01_data_preprocessing.R"))
library(glmnet)
library(caret)
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
age_sites <- readRDS(here("results/rds/age_sites.rds"))

#filter for colon samples
se <- se[,se$organ == "COL" & se$microbiome %in% c("gf", "wt")]

#filter for aging sites
se <- subsetByOverlaps(se, age_sites)

#filter for coverage/variance
se <- filter_rrbs(se, nblocks = 20, cores = 40, cpg_select = "cgi", progressbar = TRUE)

meth <- getMeth(se, type = "raw")
meth <- realize_Parallel(meth, cores = 6, nblocks = 3)
meth <- meth[!matrixStats::rowAnyNAs(meth),]

#get only the wt samples 
meth_wt <- t(meth[,se$microbiome == "wt"])
#get only the gf samples
meth_gf <- t(meth[,se$microbiome == "gf"])

age_wt <- colData(se[,se$microbiome == "wt"])$age
age_gf <- colData(se[,se$microbiome == "gf"])$age

set.seed(123)
lambda_seq <- 10^seq(2, -2, by = -.1)
#fit on wt
fit_wt <- glmnet::cv.glmnet(meth_wt, age_wt, alpha = 0)
#fit on gf
fit_gf <- glmnet::cv.glmnet(meth_gf, age_gf, alpha = 0)

#predict on wt
age_pred_wt.wt <- predict(fit_wt, newx = meth_wt, s = fit_wt$lambda.min) %>% as.double()
age_pred_gf.wt <- predict(fit_wt, newx = meth_gf, s = fit_wt$lambda.min) %>% as.double()

#predict on gf
age_pred_wt.gf <- predict(fit_gf, newx = meth_wt, s = fit_gf$lambda.min) %>% as.double()
age_pred_gf.gf <- predict(fit_gf, newx = meth_gf, s = fit_gf$lambda.min) %>% as.double()

#create dataframe with wt preds
pred_wt.wt <- data.frame(actual = age_wt, pred = age_pred_wt.wt, row.names = rownames(meth_wt))
pred_gf.wt <- data.frame(actual = age_gf, pred = age_pred_gf.wt, row.names = rownames(meth_gf))

#create dataframe with gf preds
pred_wt.gf <- data.frame(actual = age_wt, pred = age_pred_wt.gf, row.names = rownames(meth_wt))
pred_gf.gf <- data.frame(actual = age_gf, pred = age_pred_gf.gf, row.names = rownames(meth_gf))

#rbind dataframes
pred.wt <- rbind(pred_wt.wt, pred_gf.wt)
pred.gf <- rbind(pred_wt.gf, pred_gf.gf)

pred.wt <- merge(pred.wt, colData(se)[,c("microbiome", "strain", "batch", "suborgan", "sampletype")], by = "row.names") %>% as.data.frame()
pred.gf <- merge(pred.gf, colData(se)[,c("microbiome", "strain", "batch", "suborgan", "sampletype")], by = "row.names") %>% as.data.frame()

#plot predictions from wt fit
p.wt <- ggplot(data = pred.wt, mapping = aes(x = actual, y = pred, color = microbiome)) +
  geom_point(aes(shape = batch)) +
  labs(x = "Actual Age", y = "Predicted Age", title = "Ridge Regression: Actual vs Predicted Age") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + 
  geom_smooth(method = lm) +
  ggpubr::stat_regline_equation(label.y = c(30, 35), label.x = 10, show.legend = FALSE, output.type = "text")
  

ggplotly(p.wt)

#plot predictions from gf fit
p.gf <- ggplot(data = pred.gf, mapping = aes(x = actual, y = pred, color = microbiome, customdata = rownames(pred))) +
  geom_point(aes(shape = batch, color = batch)) +
  labs(x = "Actual Age", y = "Predicted Age", title = "Ridge Regression: Actual vs Predicted Age") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + 
  geom_smooth(method = lm) +
  ggpubr::stat_regline_equation(label.y = c(30, 35), label.x = 10, show.legend = FALSE, output.type = "text"); p.gf

ggplotly(p.gf)
