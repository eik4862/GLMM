source("./Code/GLMM.R")
source("./Code/NaiveBayes.R")
source("./Code/LDA.R")
source("./Code/KNN.R")
source("./Code/ANN.R")
source("./Code/RandomForest.R")
source("./Code/Logistic.R")

# Load data
load("./Data/PCA.RData")
load("./Data/model_spec.RData")
bandpow <- read.csv("./Data/bandpower_log.csv", 
                    colClasses = c("factor", "factor", "numeric", "factor", "factor", "numeric", "factor"))

# GLMM
X.glmm <- glmm.data.mat(bandpow)
glmm.model.spec <- paste0("res~", paste(colnames(X.glmm[,-(1:3)])[-c(1, 5, 18)], collapse = "+"), "+(1|subj/sess)")
glmm.cv.res <- glmm.cv(glmm.model.spec, X.glmm)

# Naive Bayes
X.naive <- naive.data.mat(bandpow)
naive.cv.res <- naive.cv(X.naive)

# LDA
X.lda <- lda.data.mat(bandpow)
lda.cv.res <- lda.cv(X.lda)

# KNN
X.knn <- knn.data.mat(bandpow, fit.pca)
knn.tune.res <- knn.tune(X.knn, pc.cand = 1:15, k.cand = 1:15)
knn.cv.res <- knn.cv(X.knn, pc = 10, k = 9)

# ANN
X.ann <- ann.data.mat(bandpow)
ann.tune.res <- ann.tune(X.ann, node.cand = 1:25, decay.cand = exp(seq(-5, 1, by = 0.25)))
ann.cv.res <- ann.cv(X.ann, node = 23, decay = exp(-0.25))

# Logistic regression
X.log <- log.data.mat(bandpow)
log.tune.res <- log.tune(X.log, link.cand = c("logit", "probit", "cauchit", "cloglog"))
log.cv.res <- log.cv(log.tune.res$formula$logit, X.log, link = "logit")

# Random Forest
X.rf <- rf.data.mat(bandpow)
rf.tune.res <- rf.tune(X.rf, mtry.cand = 1:21)
rf.cv.res <- rf.cv(X.rf, mtry = 2)

# Plots
source("./Code/PlotBase.R")

# Figure 6
# size: 6.2 * 6
plt.knn.tune <- ggplot(data = knn.tune.res$raw, aes(x = k, y = pc, fill = auc)) +
  geom_tile()+
  annotate("text", x = 9, y = 10, label = "+", size = 5) +
  labs(x = TeX("Number of neighbor ($k$)"), y = "Number of PC (d)", fill = "AUC") +
  scale_fill_viridis_c(guide = guide_colorbar(barwidth = 7, raster = T)) +
  theme_classic() +
  theme(legend.position = "bottom")

plt.ann.tune <- ggplot(data = ann.tune.res$raw, aes(x = node, y = log(decay), fill = auc)) +
  geom_tile() +
  annotate("text", x = 23, y = -0.25, label = "+", size = 5) +
  labs(x = TeX("Number of node ($n$)"), y = TeX("Decay constant ($\\log\\,d$)"), fill = "AUC") +
  scale_fill_viridis_c(guide = guide_colorbar(barwidth = 7, raster = T)) +
  theme_classic() +
  theme(legend.position = "bottom")

log.tune.res$raw$link <- factor(log.tune.res$raw$link, levels = c("logit", "probit", "cauchit", "cloglog"))
plt.log.tune <- ggplot(data = log.tune.res$raw, aes(x = link, y = auc, fill = link)) +
  geom_bar(stat = "identity") +
  labs(x = "Link function", y = "AUC") +
  annotate("text", x = log.tune.res$raw$link, y = log.tune.res$raw$auc, 
           label = round(log.tune.res$raw$auc, 3), vjust = -1) +
  coord_cartesian(ylim = c(0.745, 0.755)) +
  theme_classic() +
  theme(legend.position = "None") +
  scale_fill_npg()

plt.rf.tune <- ggplot(data = rf.tune.res$raw, aes(x = mtry, y = auc, color = "A")) +
  geom_point(data = rf.tune.res$raw[-2,], aes(x = mtry, y = auc), size = 2, shape = 3, color = "black") +
  geom_line(size = 0.8) +
  geom_point(x = 2, y = rf.tune.res$raw$auc[2], shape = 15, size = 2) +
  labs(x = TeX("Size of subset ($d$)"), y = "AUC") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg()

ggarrange(plt.knn.tune, plt.ann.tune, plt.log.tune, plt.rf.tune, labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2, heights = c(1.4, 1))

# Figure 7
# size: 9 * 6
tmp <- data.frame(stat = c(lda.cv.res$stat[1,], knn.cv.res$stat[1,], ann.cv.res$stat[1,],
                           glmm.cv.res$stat[1,], rf.cv.res$stat[1,], log.cv.res$stat[1,], naive.cv.res$stat[1,]),
           measure = rep(c("ACC", "TPR", "TNR", "F score", "MCC", "AUC"), times = 7),
           method = factor(rep(c("LDA", "KNN", "ANN", "GLMM", "RF", "LOG", "NB"), each = 6),
                           levels = c("GLMM", "NB", "LDA", "KNN", "ANN", "LOG", "RF")),
           sd = c(lda.cv.res$stat[2,], knn.cv.res$stat[2,], ann.cv.res$stat[2,], glmm.cv.res$stat[2,],
                  rf.cv.res$stat[2,], log.cv.res$stat[2,], naive.cv.res$stat[2,]))
tmp <- split(tmp, tmp$measure)

plt.acc <- ggplot(data = tmp$ACC, aes(x = method, y = stat, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = stat - sd, ymax = stat + sd), width = 0.2) +
  coord_cartesian(ylim = c(0.75, 0.9)) +
  labs(x = "", y = "ACC") +
  theme_classic() +
  theme(legend.position = "None", axis.text.x=element_blank()) +
  scale_fill_npg()

plt.auc <- ggplot(data = tmp$AUC, aes(x = method, y = stat, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = stat - sd, ymax = stat + sd), width = 0.2) +
  coord_cartesian(ylim = c(0.7, 0.9)) +
  labs(x = "", y = "AUC", fill = "Classifier") +
  theme_classic() +
  theme(legend.position = "top", axis.text.x=element_blank()) +
  guides(fill = guide_legend(byrow = T, nrow = 1)) +
  scale_fill_npg()

plt.mcc <- ggplot(data = tmp$MCC, aes(x = method, y = stat, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = stat - sd, ymax = stat + sd), width = 0.2) +
  coord_cartesian(ylim = c(0.2, 0.5)) +
  labs(x = "", y = "MCC") +
  theme_classic() +
  theme(legend.position = "None", axis.text.x=element_blank()) +
  scale_fill_npg()

plt.f <- ggplot(data = tmp$`F score`, aes(x = method, y = stat, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = stat - sd, ymax = stat + sd), width = 0.2) +
  coord_cartesian(ylim = c(0.85, 0.95)) +
  labs(x = "", y = "F score") +
  theme_classic() +
  theme(legend.position = "None", axis.text.x=element_blank()) +
  scale_fill_npg()

plt.tpr <- ggplot(data = tmp$TPR, aes(x = method, y = stat, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = stat - sd, ymax = stat + sd), width = 0.2) +
  coord_cartesian(ylim = c(0.85, 1)) +
  labs(x = "", y = "TPR") +
  theme_classic() +
  theme(legend.position = "None", axis.text.x=element_blank()) +
  scale_fill_npg()

plt.tnr <- ggplot(data = tmp$TNR, aes(x = method, y = stat, fill = method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = stat - sd, ymax = stat + sd), width = 0.2) +
  coord_cartesian(ylim = c(0.1, 0.5)) +
  labs(x = "", y = "TNR") +
  theme_classic() +
  theme(legend.position = "None", axis.text.x=element_blank()) +
  scale_fill_npg()

ggarrange(plt.auc, plt.mcc, plt.f, plt.acc, plt.tpr, plt.tnr, labels = c("A", "B", "C", "D", "E", "F"), 
          nrow = 3, ncol = 2, common.legend = T)

data.frame(tpr = c(lda.cv.res$roc$TPR, qda.cv.res$roc$TPR, knn.cv.res$roc$TPR, ann.cv.res$roc$TPR,
  glmm.cv.res$roc$TPR, rf.cv.res$roc$TPR, logi.cv.res$roc$TPR, naive.cv.res$roc$TPR),
  fpr <- rep(lda.cv.res$roc$FPR, times = 8),
  method = factor(rep(c("LDA", "QDA", "KNN", "ANN", "GLMM", "RF", "LOGI", "NB"), each = length(lda.cv.res$roc$TPR)),
                  levels = c("GLMM", "NB", "LDA", "QDA", "KNN", "ANN", "LOGI", "RF"))) %>%
  ggplot(aes(x = fpr, y = tpr, color = method)) +
  geom_line(size = 1) +
  geom_line(data = data.frame(x = c(0, 1), y = c(0, 1), method = "Ref"), 
            aes(x = x, y = y), linetype = "dashed", color = "black", size = 0.8) +
  labs(x = "FPR", y = "TPR", color = "Method") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_npg()

rm(lfp)
                    
                    
                    