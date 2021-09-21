source("./Code/ENet.R")
source("./Code/GLMM.R")

# Load data
bandpow <- read.csv("./Data/bandpower_log.csv", 
                    colClasses = c("factor", "factor", "numeric", "factor", "factor", "numeric", "factor"))

# Elastic net
X.enet <- enet.data.mat(bandpow)
enet.tune.res <- enet.tune(X.enet, alpha.cand = seq(0.01, 0.99, by = 0.01), lambda.cand = exp(seq(-1, -5, by = -0.04)), iter = 3)
fit.enet <- glmnet(as.matrix(X.enet[,-1]), X.enet$res, family = "binomial", alpha = 0.15, lambda = 0.01179594)
fit.enet$beta
colnames(X.enet)[-1]
glmm.model.spec <- paste0("res~", paste(colnames(X.enet[,-1])[-c(1, 5, 18)], collapse = "+"), "+(1|subj/sess)")
write(glmm.model.spec, "./Data/model_spec.RData")

# GLMM
X.glmm <- glmm.data.mat(bandpow)
fit.glmm <- glmer(glmm.model.spec, data = X.glmm, family = binomial,
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
round(summary(fit.glmm)$coefficients, 3)
summary(fit.glmm)

# LRT
glmm.lrt.res <- glmm.lrt(X.glmm, colnames(X.enet[,-1])[-c(1, 5, 18)])

# Plots
source("./Code/PlotBase.R")

# Figure 5
# size: 3.5 * 6
plt.enet.auc <- ggplot(enet.tune.res$raw.mean, aes(x = log(lambda), y = alpha, fill = auc)) +
  geom_raster(interpolate = T) +
  annotate("text", x = log(0.008915179), y = 0.04, label = "+", size = 5) +
  annotate("text", x = log(0.01179594), y = 0.15, label = "+", size = 5, color = "red") +
  labs(x = TeX("Penalty factor ($\\log\\,\\lambda$)"), y = TeX("Balancing factor ($\\alpha$)"), fill = "AUC") +
  scale_fill_viridis_c(guide = guide_colorbar(barwidth = 7, raster = T)) +
  theme_classic() +
  theme(legend.position = "bottom")

plt.enet.nzero <- ggplot(enet.tune.res$nzero, aes(x = log(lambda), y = alpha, fill = nzero)) +
  geom_raster(interpolate = T) +
  annotate("text", x = log(0.008915179), y = 0.04, label = "+", size = 5) +
  annotate("text", x = log(0.01179594), y = 0.15, label = "+", size = 5, color = "red") +
  labs(x = TeX("Penalty factor ($\\log\\,\\lambda$)"), y = TeX("Balancing factor ($\\alpha$)"), fill = "Nzero") +
  scale_fill_viridis_c(guide = guide_colorbar(barwidth = 7, raster = T)) +
  theme_classic() +
  theme(legend.position = "bottom")

ggarrange(plt.enet.auc, plt.enet.nzero, labels = c("A", "B"), nrow = 1, ncol = 2)

