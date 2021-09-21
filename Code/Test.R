source("./Code/Power.R")

# t.test various err dist
t.norm <- est.pow(rnorm, test = function(x, y) t.test(x, y, var.equal = TRUE),
                  m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
t.unif <- est.pow(function(n) runif(n, min = -1, max = 1), test = function(x, y) t.test(x, y, var.equal = TRUE),
                  m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
t.t <- est.pow(function(n) rt(n, df = 1), test = function(x, y) t.test(x, y, var.equal = TRUE), 
               m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
t.de <- est.pow(rlaplace, test = function(x, y) t.test(x, y, var.equal = TRUE), 
                m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))

# wilcox.test various err dist
wilcox.norm <- est.pow(rnorm, test = wilcox.test, 
                       m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
wilcox.unif <- est.pow(function(n) runif(n, min = -1, max = 1), test = wilcox.test, 
                       m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
wilcox.t <- est.pow(function(n) rt(n, df = 1), test = wilcox.test, 
                    m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
wilcox.de <- est.pow(rlaplace, test = wilcox.test, 
                     m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))

# QQ
set.seed(0)

# Sample
norm.quant <- sort(rnorm(100))
t.quant <- sort(rt(100, df = 1))
unif.quant <- sort(runif(100, min = -1, max = 1))
de.quant <- sort(rlaplace(100))
ref.quant <- sort(rnorm(100))

# Fit line
x.out <- seq(ref.quant[1], ref.quant[100], length.out = 100)
norm.qqline <- lm(y ~ x, data = data.frame(x = norm.quant[25:75], y = ref.quant[25:75]))
t.qqline <- lm(y ~ x, data = data.frame(x = t.quant[25:75], y = ref.quant[25:75]))
unif.qqline <- lm(y ~ x, data = data.frame(x = unif.quant[25:75], y = ref.quant[25:75]))
de.qqline <- lm(y ~ x, data = data.frame(x = de.quant[25:75], y = ref.quant[25:75]))


# t.test dependent sample
t.g10 <- est.pow(function(n) depend.rand(n, 10), test = function(x, y) t.test(x, y, var.equal = TRUE), 
                 m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
t.g25 <- est.pow(function(n) depend.rand(n, 4), test = function(x, y) t.test(x, y, var.equal = TRUE), 
                 m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
t.g50 <- est.pow(function(n) depend.rand(n, 2), test = function(x, y) t.test(x, y, var.equal = TRUE), 
                 m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))

# wilcox.test dependent sample
wilcox.g10 <- est.pow(function(n) depend.rand(n, 10), test = wilcox.test, m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
wilcox.g25 <- est.pow(function(n) depend.rand(n, 4), test = wilcox.test, m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))
wilcox.g50 <- est.pow(function(n) depend.rand(n, 2), test = wilcox.test, m = 100, n = 100, diff = seq(-0.5, 0.5, by = 0.01))

# SACF
set.seed(0)
acf.g10 <- as.vector(acf(depend.rand(100, 10), lag.max = 15)$acf)
acf.g25 <- as.vector(acf(depend.rand(100, 4), lag.max = 15)$acf)
acf.g50 <- as.vector(acf(depend.rand(100, 2), lag.max = 15)$acf)
acf.g100 <- as.vector(acf(depend.rand(100, 1), lag.max = 15)$acf)

# Plot
source("./Code/PlotBase.R")

# Figure 1
# size: 9 * 6
plt.t.dist <- cbind(rbind(t.unif, t.norm, t.de, t.t),
                    dist = factor(rep(c("U(-1,1)", "N(0,1)", "DE(0,1)", "t(1)"), each = nrow(t.norm)), 
                                  levels = c("U(-1,1)", "N(0,1)", "DE(0,1)", "t(1)"))) %>%
  ggplot(aes(x = diff, y = pow, color = dist)) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.1, se = F, size = 1) +
  geom_hline(yintercept = 0.05, lty = "dashed", size = 0.8) +
  labs(color = "Distribution", x = TeX("Difference in mean ($\\Delta$)"), y = TeX("Estimated power ($\\beta$)")) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_npg()

plt.wilcox.dist <- cbind(rbind(wilcox.unif, wilcox.norm, wilcox.de, wilcox.t),
                         dist = factor(rep(c("U(-1,1)", "N(0,1)", "DE(0,1)", "t(1)"), each = nrow(wilcox.norm)), 
                                       levels = c("U(-1,1)", "N(0,1)", "DE(0,1)", "t(1)"))) %>%
  ggplot(aes(x = diff, y = pow, color = dist)) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.1, se = F, size = 1) +
  geom_hline(yintercept = 0.05, lty = "dashed", size = 0.8) +
  labs(color = "Distribution", x = TeX("Difference in mean ($\\Delta$)"), y = TeX("Estimated power ($\\beta$)")) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_npg()

plt.norm.qq <- data.frame(samp = norm.quant, ref = ref.quant, x.out = x.out) %>%
  ggplot(aes(x = ref, y = samp)) +
  geom_point(size = 2, color = "#4DBBD5FF") +
  geom_segment(aes(x = ref.quant[1], xend = ref.quant[100], 
                   y = predict(norm.qqline, data.frame(x = ref.quant[1])),
                   yend = predict(norm.qqline, data.frame(x = ref.quant[100]))), size = 1, color = "#4DBBD5FF") +
  labs(x = "Theoretical", y = "Sample") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg() +
  scale_fill_npg()

plt.t.qq <- data.frame(samp = t.quant, ref = ref.quant, x.out = x.out) %>%
  ggplot(aes(x = ref, y = samp)) +
  geom_point(size = 2, color = "#3C5488FF") +
  geom_segment(aes(x = ref.quant[1], xend = ref.quant[100], 
                   y = predict(t.qqline, data.frame(x = ref.quant[1])),
                   yend = predict(t.qqline, data.frame(x = ref.quant[100]))), size = 1, color = "#3C5488FF") +
  labs(x = "Theoretical", y = "Sample") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg() +
  scale_fill_npg()

plt.unif.qq <- data.frame(samp = unif.quant, ref = ref.quant, x.out = x.out) %>%
  ggplot(aes(x = ref, y = samp)) +
  geom_point(size = 2, color = "#CD0000FF") +
  geom_segment(aes(x = ref.quant[1], xend = ref.quant[100], 
                   y = predict(unif.qqline, data.frame(x = ref.quant[1])),
                   yend = predict(unif.qqline, data.frame(x = ref.quant[100]))), size = 1, color = "#CD0000FF") +
  labs(x = "Theoretical", y = "Sample") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg() +
  scale_fill_npg()

plt.de.qq <- data.frame(samp = de.quant, ref = ref.quant, x.out = x.out) %>%
  ggplot(aes(x = ref, y = samp)) +
  geom_point(size = 2, color = "#00A087FF") +
  geom_segment(aes(x = ref.quant[1], xend = ref.quant[100], 
                   y = predict(de.qqline, data.frame(x = ref.quant[1])),
                   yend = predict(de.qqline, data.frame(x = ref.quant[100]))), size = 1, color = "#00A087FF") +
  labs(x = "Theoretical", y = "Sample") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg() +
  scale_fill_npg()

ggarrange(plt.t.dist, plt.wilcox.dist, plt.unif.qq, plt.norm.qq, plt.de.qq, plt.t.qq, 
          labels = c("A", "B", "C", "D" ,"E", "F"), common.legend = T, nrow = 3, ncol = 2)

# Figure 2
# size: 9 * 6
plt.t.dep <- cbind(rbind(t.g10, t.g25, t.g50, t.norm),
                   dist = factor(rep(c(10, 25, 50, 100), each = nrow(t.norm)), 
                                 levels = c(10, 25, 50, 100))) %>%
  ggplot(aes(x = diff, y = pow, color = dist)) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.1, se = F, size = 1) +
  geom_hline(yintercept = 0.05, lty = "dashed", size = 0.8) +
  labs(color = "# of subgroups", x = TeX("Difference in mean ($\\Delta$)"), y = TeX("Estimated power ($\\beta$)")) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_npg()

plt.wilcox.dep <- cbind(rbind(wilcox.g10, wilcox.g25, wilcox.g50, wilcox.norm),
                        dist = factor(rep(c(10, 25, 50, 100), each = nrow(wilcox.norm)), 
                                      levels = c(10, 25, 50, 100))) %>%
  ggplot(aes(x = diff, y = pow, color = dist)) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.1, se = F, size = 1) +
  geom_hline(yintercept = 0.05, lty = "dashed", size = 0.8) +
  labs(color = "# of subgroups", x = TeX("Difference in mean ($\\Delta$)"), y = TeX("Estimated power ($\\beta$)")) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_npg()

plt.acf.g10 <- data.frame(lag = 0:15, acf = acf.g10) %>%
  ggplot(aes(x = lag, y = acf)) +
  geom_segment(aes(x = 0:15, xend = 0:15, y = 0, yend = acf.g10), color = "#CD0000FF", alpha = 0.4) +
  geom_point(color = "#CD0000FF", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = c(-qnorm(0.975), ymax = qnorm(0.975)) / 10, color = "#CD0000FF", size = 0.8, lty = "dashed") +
  labs(x = "Lag", y = "SACF") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg()

plt.acf.g25 <- data.frame(lag = 0:15, acf = acf.g25) %>%
  ggplot(aes(x = lag, y = acf)) +
  geom_segment(aes(x = 0:15, xend = 0:15, y = 0, yend = acf.g25), color = "#4DBBD5FF", alpha = 0.4) +
  geom_point(color = "#4DBBD5FF", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = c(-qnorm(0.975), ymax = qnorm(0.975)) / 10, color = "#4DBBD5FF", size = 0.8, lty = "dashed") +
  labs(x = "Lag", y = "SACF") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg()

plt.acf.g50 <- data.frame(lag = 0:15, acf = acf.g50) %>%
  ggplot(aes(x = lag, y = acf)) +
  geom_segment(aes(x = 0:15, xend = 0:15, y = 0, yend = acf.g50), color = "#00A087FF", alpha = 0.4) +
  geom_point(color = "#00A087FF", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = c(-qnorm(0.975), ymax = qnorm(0.975)) / 10, color = "#00A087FF", size = 0.8, lty = "dashed") +
  labs(x = "Lag", y = "SACF") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg()

plt.acf.g100 <- data.frame(lag = 0:15, acf = acf.g100) %>%
  ggplot(aes(x = lag, y = acf)) +
  geom_segment(aes(x = 0:15, xend = 0:15, y = 0, yend = acf.g100), color = "#3C5488FF", alpha = 0.4) +
  geom_point(color = "#3C5488FF", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = c(-qnorm(0.975), ymax = qnorm(0.975)) / 10, color = "#3C5488FF", size = 0.8, lty = "dashed") +
  labs(x = "Lag", y = "SACF") +
  theme_classic() +
  theme(legend.position = "None") +
  scale_color_npg()

ggarrange(plt.t.dep, plt.wilcox.dep, plt.acf.g10, plt.acf.g25, plt.acf.g50, plt.acf.g100,
          labels = c("A", "B", "C", "D" ,"E", "F"), common.legend = T, nrow = 3, ncol = 2)
