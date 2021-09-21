source("./Code/Etc.R")

# Load data
lfp <- read.csv("./Data/LFP.csv", 
                colClasses = c("numeric", "numeric", "numeric", "factor", "factor"))
bandpow <- read.csv("./Data/bandpower.csv", 
                    colClasses = c("factor", "factor", "numeric", "factor", "factor", "numeric", "factor"))

# VST
bandpow.log <- bandpow
bandpow.log$bandpow <- log(bandpow.log$bandpow)
bandpow.log <- subset(bandpow.log, bandpow != -Inf)

# Compute hit & miss ratio
res.ratio <- ratio(bandpow)

# PCA
X.pca <- pca.data.mat(bandpow.log)
fit.pca <- pca(X.pca, scale = F)

# Plots
source("./Code/PlotBase.R")

# Figure 3
# size: 7.5 * 6
plt.lfp <- subset(lfp, chan == "AC_L") %>%
  ggplot(aes(x = time, y = signal, group = trial)) +
  geom_line(size = 0.5, alpha = 0.3) +
  geom_vline(xintercept = c(0, 1), lty = "dashed", color = "#CD0000FF", size = 0.8) +
  facet_wrap(. ~ subject, ncol = 1, strip.position = "right") +
  geom_text(data = data.frame(subject = rep("Subject 1", times = 1), 
                              text = c("Prestimulus", "Stimulus", "Poststimulus"),
                              trial = rep(1, times = 3),
                              time = c(-1, 0.5, 2), signal = rep(10, times = 3)), aes(label = text), vjust = 0.8) +
  labs(x = "Time (sec)", y = TeX("LFP signal ($\\mu V$)")) +
  theme_classic() +
  theme(legend.position = "None", strip.background = element_blank()) +
  scale_color_npg()

plt.lfp

# Figure 4
# Size: 6 * 12
plt.vst.bf <- subset(bandpow, band == "delta" & chan == "AC_L") %>%
  ggplot(aes(x = subj, y = bandpow, fill = res)) +
  geom_boxplot(outlier.size = 1.5) +
  labs(x = "Subject", y = TeX("Bandpower ($\\mu V^2/Hz$)"), fill = "Group") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_fill_npg()

plt.vst.af <- subset(bandpow.log, band == "delta" & chan == "AC_L") %>%
  ggplot(aes(x = subj, y = bandpow, fill = res)) +
  geom_boxplot(outlier.size = 1.5) +
  labs(x = "Subject", y = TeX("Bandpower ($\\mu V^2/Hz$)"), fill = "Group") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_fill_npg()

plt.kde.bf <- subset(bandpow, band == "delta" & chan == "AC_L" & subj == 2) %>%
  ggplot(aes(x = bandpow, color = res, fill = res)) +
  geom_density(size = 1) +
  geom_point(aes(x = bandpow, y = ifelse(res == "Hit", -0.16, -0.36)), size = 2, shape = "|", show.legend = F) +
  labs(x = "", y = "Estimated density", color = "Group", fill = "Group") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_npg() +
  scale_fill_npg(alpha = 0.4)

plt.kde.af <- subset(bandpow.log, band == "delta" & chan == "AC_L" & subj == 2) %>%
  ggplot(aes(x = bandpow, color = res, fill = res)) +
  geom_density(size = 1) +
  geom_point(aes(x = bandpow, y = ifelse(res == "Hit", -0.04, -0.1)), size = 2, shape = "|") +
  labs(x = "", y = "Estimated density") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg() +
  scale_fill_npg(alpha = 0.4)

bw.nrd0(subset(bandpow, band == "delta" & chan == "AC_L" & subj == 2 & res == "Hit")$bandpow)
bw.nrd0(subset(bandpow, band == "delta" & chan == "AC_L" & subj == 2 & res == "Miss")$bandpow)
bw.nrd0(subset(bandpow.log, band == "delta" & chan == "AC_L" & subj == 2 & res == "Hit")$bandpow)
bw.nrd0(subset(bandpow.log, band == "delta" & chan == "AC_L" & subj == 2 & res == "Miss")$bandpow)

plt.ratio <- ggplot(data = res.ratio, aes(x = subj, y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Subject", y = "Ratio (%)", fill = "Group") +
  theme_classic() +
  theme(legend.position = "top", strip.background = element_blank()) +
  scale_fill_npg()

plt.pc12 <- with(fit.pca, 
     data.frame(x = fit[[1]]$x[,1], y = fit[[1]]$x[,2], res = X.pca[[1]]$res) %>%
       ggplot(aes(x = x, y = y, color = res)) +
       labs(x = "PC 1", y = "PC 2") +
       geom_point(size = 2) +
       theme_classic() +
       scale_color_npg()
  )

plt.pc13 <- with(fit.pca, 
                 data.frame(x = fit[[1]]$x[,1], y = fit[[1]]$x[,3], res = X.pca[[1]]$res) %>%
                   ggplot(aes(x = x, y = y, color = res)) +
                   labs(x = "PC 1", y = "PC 3") +
                   geom_point(size = 2) +
                   theme_classic() +
                   scale_color_npg()
)

plt.pc23 <- with(fit.pca, 
                 data.frame(x = fit[[1]]$x[,2], y = fit[[1]]$x[,3], res = X.pca[[1]]$res) %>%
                   ggplot(aes(x = x, y = y, color = res)) +
                   labs(x = "PC 2", y = "PC 3") +
                   geom_point(size = 2) +
                   theme_classic() +
                   scale_color_npg()
)

ggarrange(plt.vst.bf, plt.vst.af, plt.kde.bf, plt.kde.af, plt.ratio, plt.pc12, plt.pc13, plt.pc23,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"), common.legend = T, nrow = 2, ncol = 4, legend.grob = get_legend(plt.ratio))

# 
ggplot(fit.pca$scree, aes(x = k, y = eigen, color = subj)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 5, lty = "dashed", size = 0.8) +
  labs(x = "# of PC", y = "Eigenvalue (standardized)", color = "Subject") +
  geom_point() +
  theme_classic() +
  scale_color_npg()

