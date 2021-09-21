library(class)
library(progress)
library(parallel)
library(doSNOW)
library(reshape)

source("./Code/PredBase.R")

rf.data.mat <- function(src) {
  # Remove BF chan and Mouse 1, 2, 7
  src <- subset(src, chan != "BF")
  src <- subset(src, !(subj %in% c(1, 2, 7)))
  src$subj <- droplevels(src$subj)
  src$chan <- droplevels(src$chan)
  
  # Reform dataframe
  data <- cast(src, subj + sess + trial + res ~ chan + band, value = "bandpow")
  
  # Remove NAs
  data <- na.omit(data)
  
  # Drop trial
  data <- data[,-c(2, 3)]
  
  return(data)
}

rf.cv <- function(data.mat, mtry, fold = 10, iter = 100, prob.cut = seq(0, 1, by = 0.01), x.out = seq(0, 1, by = 0.01), seed = T, parallel = T, verbose = T) {
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
    
    blk.size <- round(nrow(data.mat) / fold)
    prog <- function(n) if (verbose) setTxtProgressBar(txtProgressBar(max = iter, style = 3), n)
    
    if (verbose) message("Performing CV")
    
    iter.res <- foreach(i = seq_len(iter), .combine = list, .multicombine = T, 
                        .packages = "randomForest", .export = c("acc", "tpr.hit", "tpr.miss", "fscore", "mcc", "auc"),
                        .options.snow = list(progress = prog)) %dopar% {  # For each iteration
      stats.tmp <- data.frame(acc = vector(length = fold), tpr.hit = vector(length = fold), tpr.miss = vector(length = fold),
                              fscore = vector(length = fold), mcc = vector(length = fold))
      roc.tpr.tmp <- vector(length = length(prob.cut))
      roc.fpr.tmp <- vector(length = length(prob.cut))
      roc.tmp <- matrix(nrow = fold, ncol = length(x.out))
      
      # Permute samples
      if (seed) set.seed(seed.list[i])
      idx.perm <- sample(nrow(data.mat))
      
      for (j in seq_len(fold)) {  # For each fold
        # Determine hold-out index
        if (j == fold) {
          hold.idx <- idx.perm[(blk.size * (j - 1) + 1):length(idx.perm)]
        } else {
          hold.idx <- idx.perm[(blk.size * (j - 1) + 1):(blk.size * j)]
        }
        
        # Train
        fit <- randomForest(res ~ ., data = data.mat[-hold.idx,], mtry = mtry)
        
        # Test
        pred <- predict(fit, data.mat[hold.idx,], type = "prob")
        confu.mat <- table(pred = factor(ifelse(pred[,1] > 0.5, "Hit", "Miss"), levels = c("Hit", "Miss")),
                           real = factor(data.mat[hold.idx,]$res, levels = c("Hit", "Miss")))
        
        # Compute stats
        stats.tmp$acc[j] <- acc(confu.mat)
        stats.tmp$tpr.hit[j] <- tpr.hit(confu.mat)
        stats.tmp$tpr.miss[j] <- tpr.miss(confu.mat)
        stats.tmp$fscore[j] <- fscore(confu.mat)
        stats.tmp$mcc[j] <- mcc(confu.mat)
        
        # Compute ROC
        for (k in seq_along(prob.cut)) {
          confu.mat <- table(pred = factor(ifelse(pred[,1] > prob.cut[k], "Hit", "Miss"), levels = c("Hit", "Miss")),
                             real = factor(data.mat[hold.idx,]$res, levels = c("Hit", "Miss")))
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
    
    blk.size <- round(nrow(data.mat) / fold)
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
      idx.perm <- sample(nrow(data.mat))
      
      for (j in seq_len(fold)) {  # For each fold
        # Determine hold-out index
        if (j == fold) {
          hold.idx <- idx.perm[(blk.size * (j - 1) + 1):length(idx.perm)]
        } else {
          hold.idx <- idx.perm[(blk.size * (j - 1) + 1):(blk.size * j)]
        }
        
        # Train
        fit <- randomForest(res ~ ., data = data.mat[-hold.idx,], mtry = mtry)
        
        # Test
        pred <- predict(fit, data.mat[hold.idx,], type = "prob")
        confu.mat <- table(pred = factor(ifelse(pred[,1] > 0.5, "Hit", "Miss"), levels = c("Hit", "Miss")),
                           real = factor(data.mat[hold.idx,]$res, levels = c("Hit", "Miss")))
        
        # Compute stats
        stats.tmp$acc[j] <- acc(confu.mat)
        stats.tmp$tpr.hit[j] <- tpr.hit(confu.mat)
        stats.tmp$tpr.miss[j] <- tpr.miss(confu.mat)
        stats.tmp$fscore[j] <- fscore(confu.mat)
        stats.tmp$mcc[j] <- mcc(confu.mat)
        
        # Compute ROC
        for (k in seq_along(prob.cut)) {
          confu.mat <- table(pred = factor(ifelse(pred[,1] > prob.cut[k], "Hit", "Miss"), levels = c("Hit", "Miss")),
                             real = factor(data.mat[hold.idx,]$res, levels = c("Hit", "Miss")))
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
    cat("\n\tCross Validation Procedure for Random Forest\n\n")
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

rf.tune <- function(data.mat, mtry.cand, fold = 10, iter = 5, seed = T, parallel = T, verbose = T) {
  # Set seed if needed
  if (seed) set.seed(0)
  
  # Setting up
  if (verbose) message("Setting up for tuning")
  
  start <- Sys.time()
  cand <- data.frame(mtry = mtry.cand)
  
  if (parallel) { # Do parallel
    if (seed) {
      seed.list <- round(1000 * rnorm(nrow(cand)))
    }
    
    cl <- makeCluster(detectCores())
    registerDoSNOW(cl)
    
    prog <- function(n) if (verbose) setTxtProgressBar(txtProgressBar(max = nrow(cand), style = 3), n)
    
    # Tune
    if (verbose) message("Performing CV")
    iter.res <- foreach (i = seq_len(nrow(cand)), .combine = rbind, .packages = "randomForest",
                         .export = c("acc", "tpr.hit", "tpr.miss", "fscore", "mcc", "auc", "rf.cv"),
                         .options.snow = list(progress = prog)) %dopar% {  # For all combination
      if (seed) set.seed(seed.list[i])
     
      rf.cv(data.mat, cand$mtry[i], fold = fold, iter = iter, seed = F, parallel = F, verbose = F)$stats[1,]
    }
    
    if (verbose) cat("\n")
    stopCluster(cl)
    
    # join
    cand <- cbind(cand, iter.res)
  } else {  # Do sequential
    cand$acc <- vector(length = nrow(cand))
    cand$tpr.hit <- vector(length = nrow(cand))
    cand$tpr.miss <- vector(length = nrow(cand))
    cand$fscore <- vector(length = nrow(cand))
    cand$mcc <- vector(length = nrow(cand))
    cand$auc <- vector(length = nrow(cand))
    
    if (verbose) {
      message("Performing CV")
      pb <- progress_bar$new(total = nrow(cand), clear = F)
    }
    
    for (i in seq_len(nrow(cand))) {  # For all combination
      # Tune
      if (verbose) pb$tick()
      cand[i,2:7] <- rf.cv(data.mat, cand$mtry[i], fold = fold, iter = iter, seed = F, parallel = F, verbose = F)$stats[1,]
    }
  }
  
  # Select best one
  if (verbose) message("Tuning hyperparameter")
  
  max.idx <- apply(cand[,2:7], MARGIN = 2, which.max)
  max.val <- c(cand$acc[max.idx[1]], cand$tpr.hit[max.idx[2]], cand$tpr.miss[max.idx[3]], 
               cand$fscore[max.idx[4]], cand$mcc[max.idx[5]], cand$auc[max.idx[6]])
  sel <- t(data.frame(mtry = cand$mtry[max.idx], max = max.val))
  colnames(sel) <- c("ACC", "TPR(Hit)", "TPR(Miss)", "F score", "MCC", "AUC")
  rownames(sel) <- c("mtry", "max")
  elapsed <- difftime(Sys.time(), start, units = "secs")[[1]]
  
  # Report
  if (verbose) {
    cat("\n\tHyperparameter Tuning Procedure for Random Forest\n\n")
    cat(sprintf("  @data    : %s\n", deparse(substitute(data.mat))))
    cat(sprintf("  @fold    : %d\n", fold))
    cat(sprintf("  @iter    : %d\n", iter))
    cat(sprintf("  @parallel: %s\n", as.character(parallel)))
    cat(sprintf("  @seed    : %s\n", as.character(seed)))
    cat(sprintf("  @elapsed : %.4f (sec)\n", elapsed))
    
    cat("\n\tTuning Results\n\n")
    print(round(sel, 4))
  }
  
  return(invisible(list(select = sel, raw = cand)))
}