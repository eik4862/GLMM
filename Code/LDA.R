library(MASS)
library(progress)
library(parallel)
library(doSNOW)
library(reshape)

source("./Code/PredBase.R")

lda.data.mat <- function(src) {
  # Reform dataframe
  data <- cast(src, subj + sess + trial + res ~ chan + band, value = "bandpow")
  
  # Drop sess & trial
  data <- data[,-(2:3)]
  
  # Split for each mouse
  data <- split(data, data$subj, drop = TRUE)
  
  # Remove NAs
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

lda.cv <- function(data.mat, fold = 10, iter = 100, prob.cut = seq(0, 1, by = 0.01), x.out = seq(0, 1, by = 0.01),
                   seed = T, parallel = T, verbose = T) {
  # Set seed if needed
  if (seed) set.seed(0)
  
  start <- Sys.time()
  
  if (parallel) { # Do parallel
    # Settin up
    if (verbose) message("Setting up for CV")
    
    if (seed) {
      seed.list <- round(1000 * rnorm(iter))
    }
    
    cl <- makeCluster(detectCores())
    registerDoSNOW(cl)
    
    blk.size <- sapply(data.mat, function(x) round(nrow(x) / fold))
    prog <- function(n) if (verbose) setTxtProgressBar(txtProgressBar(max = iter, style = 3), n)
    
    if (verbose) message("Performing CV")
    
    iter.res <- foreach(i = seq_len(iter), .combine = list, .multicombine = T, 
                        .packages = "MASS", .export = c("acc", "tpr.hit", "tpr.miss", "fscore", "mcc", "auc"),
                        .options.snow = list(progress = prog)) %dopar% {  # For each iteration
      stats.tmp <- data.frame(acc = vector(length = fold), tpr.hit = vector(length = fold), tpr.miss = vector(length = fold),
                              fscore = vector(length = fold), mcc = vector(length = fold))
      roc.tpr.tmp <- vector(length = length(prob.cut))
      roc.fpr.tmp <- vector(length = length(prob.cut))
      roc.tmp <- matrix(nrow = fold, ncol = length(x.out))
      
      # Permute samples
      if (seed) set.seed(seed.list[i])
      idx.perm <- sapply(data.mat, function(x) sample(nrow(x)))
      
      for (j in seq_len(fold)) {  # For each fold
        # Determine hold-out index
        if (j == fold) {
          hold.idx <- lapply(seq_along(data.mat), function(k) idx.perm[[k]][(blk.size[k] * (j - 1) + 1):length(idx.perm[[k]])])
        } else {
          hold.idx <- lapply(seq_along(data.mat), function(k) idx.perm[[k]][(blk.size[k] * (j - 1) + 1):(blk.size[k] * j)])
        }
        
        # Train
        pred <- lapply(seq_along(data.mat), function(k) {
          predict(lda(res ~ ., data = data.mat[[k]][-hold.idx[[k]],]), data.mat[[k]][hold.idx[[k]],])
        })
        
        # Test
        confu.mat <- lapply(seq_along(data.mat), function(k) {
          table(pred = factor(pred[[k]]$class, level = c("Hit", "Miss")), 
                real = factor(data.mat[[k]]$res[hold.idx[[k]]], levels = c("Hit", "Miss")))
        })
        confu.mat <- Reduce("+", confu.mat)
        
        # Compute stats
        stats.tmp$acc[j] <- acc(confu.mat)
        stats.tmp$tpr.hit[j] <- tpr.hit(confu.mat)
        stats.tmp$tpr.miss[j] <- tpr.miss(confu.mat)
        stats.tmp$fscore[j] <- fscore(confu.mat)
        stats.tmp$mcc[j] <- mcc(confu.mat)
        
        # Compute ROC
        for (k in seq_along(prob.cut)) {
          confu.mat <- lapply(seq_along(data.mat), function(l) {
            table(pred = factor(ifelse(pred[[l]]$posterior[,1] > prob.cut[k], "Hit", "Miss"), levels = c("Hit", "Miss")),
                  real = factor(data.mat[[l]]$res[hold.idx[[l]]], levels = c("Hit", "Miss")))
          })
          confu.mat <- Reduce("+", confu.mat)
          roc.tpr.tmp[k] <- confu.mat[1,1] / sum(confu.mat[,1])
          roc.fpr.tmp[k] <- confu.mat[1,2] / sum(confu.mat[,2])
        }
        
        options(warn = -1)
        roc.tmp[j,] <- approx(roc.fpr.tmp, roc.tpr.tmp, x.out, method = "linear")$y
        options(warn = 0)
      }
      
      # Mean over folds
      stats.raw <- colMeans(stats.tmp, na.rm = T)
      roc.raw <- colMeans(roc.tmp)
      stats.raw[6] <- auc(roc.raw)
      
      return(list(stats.raw = stats.raw, roc.raw = roc.raw))
    }
    
    if (verbose) cat("\n")
    stopCluster(cl)
    
    # Join
    stats.raw <- t(sapply(iter.res, function(x) x$stats.raw))
    roc.raw <- t(sapply(iter.res, function(x) x$roc.raw))
    colnames(stats.raw)[6] <- "auc"
  } else {  # Do sequential
    # Setting up
    if (verbose) {
      message("Setting up for CV")
      pb <- progress_bar$new(total = iter, clear = F)
    }
    
    blk.size <- sapply(data.mat, function(x) round(nrow(x) / fold))
    stats.raw <- data.frame(acc = vector(length = iter), tpr.hit = vector(length = iter), tpr.miss = vector(length = iter),
                            fscore = vector(length = iter), mcc = vector(length = iter), auc = vector(length = iter))
    roc.raw <- matrix(nrow = iter, ncol = length(x.out))
    stats.tmp <- data.frame(acc = vector(length = fold), tpr.hit = vector(length = fold), tpr.miss = vector(length = fold),
                            fscore = vector(length = fold), mcc = vector(length = fold))
    roc.tpr.tmp <- vector(length = length(prob.cut))
    roc.fpr.tmp <- vector(length = length(prob.cut))
    roc.tmp <- matrix(nrow = fold, ncol = length(x.out))
    
    if (verbose) message("Performing CV")
    
    for (i in seq_len(iter)) {  # For each iteration
      # Permute samples
      if (verbose) pb$tick()
      idx.perm <- sapply(data.mat, function(x) sample(nrow(x)))
      
      for (j in seq_len(fold)) {  # For each fold
        # Determine hold-out index
        if (j == fold) {
          hold.idx <- lapply(seq_along(data.mat), function(k) idx.perm[[k]][(blk.size[k] * (j - 1) + 1):length(idx.perm[[k]])])
        } else {
          hold.idx <- lapply(seq_along(data.mat), function(k) idx.perm[[k]][(blk.size[k] * (j - 1) + 1):(blk.size[k] * j)])
        }
        
        # Train
        pred <- lapply(seq_along(data.mat), function(k) {
          predict(lda(res ~ ., data = data.mat[[k]][-hold.idx[[k]],]), data.mat[[k]][hold.idx[[k]],])
        })
        
        # Test
        confu.mat <- lapply(seq_along(data.mat), function(k) {
          table(pred = factor(pred[[k]]$class, level = c("Hit", "Miss")), 
                real = factor(data.mat[[k]]$res[hold.idx[[k]]], levels = c("Hit", "Miss")))
        })
        confu.mat <- Reduce("+", confu.mat)
        
        # Compute stats
        stats.tmp$acc[j] <- acc(confu.mat)
        stats.tmp$tpr.hit[j] <- tpr.hit(confu.mat)
        stats.tmp$tpr.miss[j] <- tpr.miss(confu.mat)
        stats.tmp$fscore[j] <- fscore(confu.mat)
        stats.tmp$mcc[j] <- mcc(confu.mat)
        
        # Compute ROC
        for (k in seq_along(prob.cut)) {
          confu.mat <- lapply(seq_along(data.mat), function(l) {
            table(pred = factor(ifelse(pred[[l]]$posterior[,1] > prob.cut[k], "Hit", "Miss"), levels = c("Hit", "Miss")),
                  real = factor(data.mat[[l]]$res[hold.idx[[l]]], levels = c("Hit", "Miss")))
          })
          confu.mat <- Reduce("+", confu.mat)
          roc.tpr.tmp[k] <- confu.mat[1,1] / sum(confu.mat[,1])
          roc.fpr.tmp[k] <- confu.mat[1,2] / sum(confu.mat[,2])
        }
        
        options(warn = -1)
        roc.tmp[j,] <- approx(roc.fpr.tmp, roc.tpr.tmp, x.out, method = "linear")$y
        options(warn = 0)
      }
      
      # Mean over folds
      stats.raw[i,-6] <- colMeans(stats.tmp, na.rm = T)
      roc.raw[i,] <- colMeans(roc.tmp)
      stats.raw$auc[i] <- auc(roc.raw[i,])
    }
  }
  
  # Compute summary stats
  if (verbose) message("Computing summary statistics")
  
  stats <- rbind(colMeans(stats.raw, na.rm = T), apply(stats.raw, MARGIN = 2, FUN = sd, na.rm = T))
  roc.mean <- colMeans(roc.raw)
  roc.sd <- apply(roc.raw, MARGIN = 2, FUN = sd)
  roc <- data.frame(fpr = x.out, tpr = roc.mean, lwr = roc.mean - roc.sd, upr = roc.mean + roc.sd)
  colnames(stats) <- c("ACC", "TPR(Hit)", "TPR(Miss)", "F score", "MCC", "AUC")
  rownames(stats) <- c("mean", "sd")
  colnames(roc) <- c("FPR", "TPR", "lower", "upper")
  elapsed <- difftime(Sys.time(), start, units = "secs")[[1]]
  
  # Report
  if (verbose) {
    cat("\n\tCross Validation Procedure for LDA\n\n")
    cat(sprintf("  @data    : %s\n", deparse(substitute(data.mat))))
    cat(sprintf("  @fold    : %d\n", fold))
    cat(sprintf("  @iter    : %d\n", iter))
    cat(sprintf("  @parallel: %s\n", as.character(parallel)))
    cat(sprintf("  @seed    : %s\n", as.character(seed)))
    cat(sprintf("  @elapsed : %.4f (sec)\n", elapsed))
    
    cat("\n\tSummary Statistics\n\n")
    print(round(stats, 4))
    
  plot(roc$FPR, roc$TPR, main = "Estimated ROC", xlab = "FPR", ylab = "TPR", sub = sprintf("AUC: %4f", stats[1,6]), type = "l", col = 2)
    lines(roc$FPR, roc$lower, lty = 2)
    lines(roc$FPR, roc$upper, lty = 2)
  }
  
  return(invisible(list(stat = stats, roc = roc, stats.raw = stats.raw, roc.raw = roc.raw)))
}
