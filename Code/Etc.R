library(progress)
library(reshape)

cutoff <- function(data, ratio = 0.03, verbose = T) {
  # Setting up
  if (verbose) message("Setting up")
  start <- Sys.time()
  cnt.bf <- nrow(data)
  
  # Split
  data$ref <- paste(data$mouse, data$chan, data$band, data$res, sep = "_")
  data.split <- split(data, data$ref, drop = T)
  
  # Cut off
  if (verbose) message("Cutting off")
  data.split <- lapply(data.split, function(x) x[order(x$bandpow),][max(1, round(ratio * nrow(x))):round((1 - ratio) * nrow(x)),])
  
  # Merge
  if (verbose) message("Merging")
  data <- Reduce(rbind, data.split)
  data <- data[,-8]
  
  cnt.af <- nrow(data)
  elapsed <- difftime(Sys.time(), start, units = "secs")[[1]]
  
  # Report
  if (verbose) {
    cat("\n\tCutting off extreme bandpowers\n\n")
    cat(sprintf("  @ratio  : %.4f\n", ratio))
    cat(sprintf("  @removed: %d\n", cnt.bf - cnt.af))
    cat(sprintf("  @elapsed: %.4f (sec)\n", elapsed))
  }
  
  return(invisible(data))
}

ratio <- function(data) {
  # Find pivot
  data <- subset(bandpow, chan == "AC_R" & band == "beta")
  
  # Compute ratio
  data$ref.hit <- data$res == "Hit"
  data$ref.miss <- data$res == "Miss"
  data$ref.tot <- rep(1, times = nrow(data))
  cnt <- aggregate(data[,8:10], by = list(data$subj), FUN = sum)
  cnt <- data.frame(subj = rep(cnt[,1], times = 2), cnt = c(cnt[,2], cnt[,3]),
                    ratio = 100 * c(cnt[,2], cnt[,3]) / rep(cnt[,4], times = 2), 
                    type = rep(c("Hit", "Miss"), each = nrow(cnt)))
  
  return(cnt)
}

pca.data.mat <- function(src, verbose = T) {
  # Reform dataframe
  if (verbose) message("Reforming data frame")
  data <- cast(src, subj + sess + trial + res ~ chan + band, value = "bandpow")
  
  # Drop sess & trial
  if (verbose) message("Dropping redundant info")
  data <- data[,-(2:3)]
  
  # Split for each mouse
  if (verbose) message("Splitting")
  data <- split(data, data$subj, drop = TRUE)
  
  # Remove NAs
  if (verbose) message("Removing NAs")
  data[[1]] <- na.omit(data[[1]][,-c(1, 3:7)])
  data[[2]] <- na.omit(data[[2]][,-c(1, 18:22)])
  data[[3]] <- na.omit(data[[3]][,-c(1, 13:17)])
  data[[4]] <- na.omit(data[[4]][,-c(1, 13:17)])
  data[[5]] <- na.omit(data[[5]][,-c(1, 13:17)])
  data[[6]] <- na.omit(data[[6]][,-1])
  data[[7]] <- na.omit(data[[7]][,-c(1, 3:7, 13:17)])
  data[[8]] <- na.omit(data[[8]][,-1])
  
  return(data)
}

pca <- function(data.mat, thres = 0.9, center = T, scale = F, verbose = T) {
  # Setting up
  if (verbose) message("Setting up")
  start <- Sys.time()
  
  # Perform PCA
  if (verbose) message("Performing PCA")
  fit <- lapply(data.mat, function(x) prcomp(x[,-1], center = center, scale = scale))
  
  # Determining optimal # of PCs
  if (verbose) message("Determining optimal # of PCs")
  # Scree plot
  eigen <- lapply(fit, function(x) x$sdev ** 2)
  scree <- data.frame(eigen = unlist(lapply(eigen, function(x) (x - min(x)) / (max(x) - min(x)))), 
                      k = c(rep(1:20, times = 5), 1:25, 1:15, 1:25),
                      subj = c(rep(1:5, each = 20), rep(6, times = 25), rep(7, times = 15), rep(8, times = 25)))
  scree$subj <- factor(scree$subj)
  
  # Kaiser
  eigen.mean <- data.frame(mean = sapply(eigen, mean))
  kaiser.pc <- sapply(seq_along(eigen.mean), function(i) sum(eigen[[i]] > eigen.mean[i]))
  
  # Cumulative variance
  cum.var <- lapply(fit, function(x) summary(x)$importance[3,])
  cum.var.pc <- sapply(fit, function(x) 1 + sum(summary(x)$importance[3,] < thres))
  cum.var <- sapply(seq_along(cum.var), function(i) cum.var[[i]][cum.var.pc[i]]) * 100
  
  pc <- t(data.frame(kaiser = kaiser.pc, mean = eigen.mean, cum.var = cum.var.pc, explained = cum.var))
  colnames(pc) <- paste("subj", 1:8)
  
  elapsed <- Sys.time() - start
  
  # Report
  if (verbose) {
    cat("\n\tPrincipal Component Analysis\n\n")
    cat(sprintf("  @data     : %s\n", deparse(substitute(data.mat))))
    cat(sprintf("  @center   : %s\n", as.character(center)))
    cat(sprintf("  @scale    : %s\n", as.character(scale)))
    cat(sprintf("  @threshold: %.4f\n", thres))
    cat(sprintf("  @elapsed  : %.4f (sec)\n", elapsed))
    
    cat("\n\tOptimal # of PCs\n\n")
    print(round(pc, 4))
  }
  
  return(invisible(list(pc = pc, fit = fit, scree = scree)))
}