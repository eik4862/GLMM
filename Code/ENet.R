library(glmnet)
library(parallel)
library(progress)
library(doSNOW)
library(reshape)

source("./Code/PredBase.R")

enet.data.mat <- function(src) {
  # Remove BF chan and Mouse 1, 2, 7
  src <- subset(src, chan != "BF")
  src <- subset(src, !(subj %in% c(1, 2, 7)))
  src$subj <- droplevels(src$subj)
  src$chan <- droplevels(src$chan)
  
  # Reform dataframe
  data <- cast(src, subj + sess + trial + res ~ chan + band, value = "bandpow")
  
  # Remove NAs
  data <- na.omit(data)
  
  # Drop mouse & sess & trial
  data <- data[,-(1:3)]

  return(data)
}

enet.cv <- function(data.mat, alpha, lambda, fold = 10, iter = 100, prob.cut = seq(0, 1, by = 0.01), x.out = seq(0, 1, by = 0.01)) {
  # Setting up
  blk.size <- round(nrow(data.mat) / fold)
  stats.raw <- data.frame(acc = vector(length = iter * length(lambda)), tpr.hit = vector(length = iter * length(lambda)),
                          tpr.miss = vector(length = iter * length(lambda)), fscore = vector(length = iter * length(lambda)),
                          mcc = vector(length = iter * length(lambda)), auc = vector(length = iter * length(lambda)),
                          lambda = factor(rep(seq_along(lambda), times = iter)))
  roc.raw <- data.frame(y = vector(length = length(x.out) * length(lambda)),
                        lambda = factor(rep(seq_along(lambda), each = length(x.out))))
  stats.tmp <- data.frame(acc = vector(length = fold * length(lambda)), tpr.hit = vector(length = fold * length(lambda)),
                          tpr.miss = vector(length = fold * length(lambda)), fscore = vector(length = fold * length(lambda)), 
                          mcc = vector(length = fold * length(lambda)), lambda = factor(rep(seq_along(lambda), times = fold)))
  roc.tpr.tmp <- matrix(ncol = length(prob.cut), nrow = length(lambda))
  roc.fpr.tmp <- matrix(ncol = length(prob.cut), nrow = length(lambda))
  roc.tmp <- data.frame(x = rep(x.out, times = length(lambda) * fold), y = vector(length = length(x.out) * length(lambda) * fold),
                        lambda = factor(rep(rep(seq_along(lambda), each = length(x.out)), times = fold)))
  
  for (i in seq_len(iter)) {  # For each iteration
    # Permute samples
    idx.perm <- sample(nrow(data.mat))
    
    for (j in seq_len(fold)) {  # For each fold
      # Determine hold-out index
      if (j == fold) {
        hold.idx <- idx.perm[(blk.size * (j - 1) + 1):length(idx.perm)]
      } else {
        hold.idx <- idx.perm[(blk.size * (j - 1) + 1):(blk.size * j)]
      }
      
      # Train
      fit <- glmnet(as.matrix(data.mat[-hold.idx,-1]), data.mat[-hold.idx,1], family = "binomial", alpha = alpha, lambda = lambda)
      
      # Test
      pred <- predict(fit, as.matrix(data.mat[hold.idx,-1]), type = "response")
      pred <- split(pred, rep(1:ncol(pred), each = nrow(pred)))
      confu.mat <- lapply(pred, function(x) table(pred = factor(ifelse(x < 0.5, "Hit", "Miss"), levels = c("Hit", "Miss")),
                                                  real = factor(data.mat[hold.idx,]$res, levels = c("Hit", "Miss"))))
      
      # Compute stats
      stats.tmp$acc[((j - 1) * length(lambda) + 1):(j * length(lambda))] <- sapply(confu.mat, acc)
      stats.tmp$tpr.hit[((j - 1) * length(lambda) + 1):(j * length(lambda))] <- sapply(confu.mat, tpr.hit)
      stats.tmp$tpr.miss[((j - 1) * length(lambda) + 1):(j * length(lambda))] <- sapply(confu.mat, tpr.miss)
      stats.tmp$fscore[((j - 1) * length(lambda) + 1):(j * length(lambda))] <- sapply(confu.mat, fscore)
      stats.tmp$mcc[((j - 1) * length(lambda) + 1):(j * length(lambda))] <- sapply(confu.mat, mcc)

      # Compute ROC
      for (k in seq_along(prob.cut)) {
        confu.mat <- lapply(pred, function(x) table(pred = factor(ifelse(x < prob.cut[k], "Hit", "Miss"), levels = c("Hit", "Miss")),
                                                    real = factor(data.mat[hold.idx,]$res, levels = c("Hit", "Miss"))))
        roc.tpr.tmp[,k] <- sapply(confu.mat, function(x) x[1,1] / sum(x[,1]))
        roc.fpr.tmp[,k] <- sapply(confu.mat, function(x) x[1,2] / sum(x[,2]))
      }
      
      options(warn = -1)
      roc.tmp$y[((j - 1) * length(x.out) * length(lambda) + 1):(j * length(x.out) * length(lambda))] <- as.vector(sapply(seq_along(lambda), 
                                                                                                                         function(k) approx(roc.fpr.tmp[k,], roc.tpr.tmp[k,], x.out, method = "linear")$y))
      options(warn = 0)
    }
    
    # Mean over folds
    stats.raw[((i - 1) * length(lambda) + 1):(i * length(lambda)),-c(6, 7)] <- aggregate(.~lambda, data = stats.tmp, FUN = mean, na.rm = T, na.action = NULL)[,-1]
    roc.raw$y <- as.vector(sapply(split(roc.tmp, roc.tmp$lambda), function(data) aggregate(.~x, data = data, FUN = mean, na.rm = T, drop = F)$y))
    stats.raw$auc[((i - 1) * length(lambda) + 1):(i * length(lambda))] <- sapply(split(roc.raw, roc.raw$lambda), function(data) auc(data$y))
  }
  
  # Compute summary stats
  stats.mean <- aggregate(.~lambda, data = stats.raw, FUN = mean, na.rm = T, na.action = NULL)
  stats.sd <- aggregate(.~lambda, data = stats.raw, FUN = sd, na.rm = T, na.action = NULL)
  
  # Compuse sparsity for full model
  nzero <- glmnet(as.matrix(data.mat[,-1]), data.mat$res, family = "binomial", alpha = alpha, lambda = lambda)$df
  
  return(invisible(list(mean = stats.mean, sd = stats.sd, nzero = nzero)))
}

enet.tune <- function(data.mat, alpha.cand, lambda.cand, fold = 10, iter = 5, seed = T, parallel = T, verbose = T) {
  # Set seed if needed
  if (seed) set.seed(0)
  
  # Setting up
  if (verbose) message("Setting up for tuning")
  
  start <- Sys.time()
  cand.mean <- data.frame(alpha = rep(alpha.cand, each = length(lambda.cand)), lambda = rep(lambda.cand, times = length(alpha.cand)))
  cand.sd <- data.frame(alpha = rep(alpha.cand, each = length(lambda.cand)), lambda = rep(lambda.cand, times = length(alpha.cand)))
  cand.nzero <- data.frame(alpha = rep(alpha.cand, each = length(lambda.cand)), lambda = rep(lambda.cand, times = length(alpha.cand)))
  
  if (parallel) { # Do parallel
    if (seed) {
      seed.list <- round(1000 * rnorm(length(alpha.cand)))
    }
    
    cl <- makeCluster(detectCores())
    registerDoSNOW(cl)
    
    prog <- function(n) if (verbose) setTxtProgressBar(txtProgressBar(max = length(alpha.cand), style = 3), n)
    
    # Tune
    if (verbose) message("Performing CV")
    iter.res <- foreach (i = seq_along(alpha.cand), .combine = list, .multicombine = T, .packages = "glmnet",
                         .export = c("acc", "tpr.hit", "tpr.miss", "fscore", "mcc", "auc", "enet.cv"),
                         .options.snow = list(progress = prog)) %dopar% {  # For all combination
                           if (seed) set.seed(seed.list[i])
                           
                           return(enet.cv(data.mat, alpha.cand[i], lambda.cand, fold = fold, iter = iter))
                         }
    
    if (verbose) cat("\n")
    stopCluster(cl)
    
    # join
    cand.mean <- cbind(cand.mean, Reduce(rbind, lapply(iter.res, function(x) x$mean))[,-1])
    cand.sd <- cbind(cand.sd, Reduce(rbind, lapply(iter.res, function(x) x$sd))[,-1])
    cand.nzero$nzero <- as.vector(sapply(iter.res, function(x) x$nzero))
  } else {  # Do sequential
    cand.mean$acc <- vector(length = nrow(cand.mean))
    cand.mean$tpr.hit <- vector(length = nrow(cand.mean))
    cand.mean$tpr.miss <- vector(length = nrow(cand.mean))
    cand.mean$fscore <- vector(length = nrow(cand.mean))
    cand.mean$mcc <- vector(length = nrow(cand.mean))
    cand.mean$auc <- vector(length = nrow(cand.mean))
    cand.sd$acc <- vector(length = nrow(cand.sd))
    cand.sd$tpr.hit <- vector(length = nrow(cand.sd))
    cand.sd$tpr.miss <- vector(length = nrow(cand.sd))
    cand.sd$fscore <- vector(length = nrow(cand.sd))
    cand.sd$mcc <- vector(length = nrow(cand.sd))
    cand.sd$auc <- vector(length = nrow(cand.sd))
    cand.nzero$nzero <- vector(length = nrow(cand.nzero))

    if (verbose) {
      message("Performing CV")
      pb <- progress_bar$new(total = nrow(alpha.cand), clear = F)
    }
    
    for (i in seq_along(alpha.cand)) {  # For all combination
      # Tune
      if (verbose) pb$tick()
      cv.res <- enet.cv(data.mat, alpha.cand[i], lambda.cand, fold = fold, iter = iter)
      cand.mean[((i - 1) * length(lambda) + 1):(i * length(lambda)),3:9] <- cv.res$mean[,-1]
      cand.sd[((i - 1) * length(lambda) + 1):(i * length(lambda)),3:9] <- cv.res$sd[,-1]
      cand.nzero$nzero[((i - 1) * length(lambda) + 1):(i * length(lambda))] <- cv.res$nzero
    }
  }
  
  # Select best one
  if (verbose) message("Tuning hyperparameter")
  
  max.idx <- apply(cand.mean[,3:8], MARGIN = 2, which.max)
  max.val <- with(cand.mean, c(acc[max.idx[1]], tpr.hit[max.idx[2]], tpr.miss[max.idx[3]], 
                               fscore[max.idx[4]], mcc[max.idx[5]], auc[max.idx[6]]))
  lwr <- max.val - with(cand.sd, c(acc[max.idx[1]], tpr.hit[max.idx[2]], tpr.miss[max.idx[3]], 
                                   fscore[max.idx[4]], mcc[max.idx[5]], auc[max.idx[6]]))
  sd1.idx <- sapply(1:6, function(i) {
    cand.msk <- cand.mean[,2 + i] >= lwr[i]
    
    return((1:nrow(cand.mean))[cand.msk][which.min(cand.nzero$nzero[cand.msk])])
  })
  sd1.val <- with(cand.mean, c(acc[sd1.idx[1]], tpr.hit[sd1.idx[2]], tpr.miss[sd1.idx[3]], 
                               fscore[sd1.idx[4]], mcc[sd1.idx[5]], auc[sd1.idx[6]]))
  max.sel <- t(data.frame(alpha = cand.mean$alpha[max.idx], lambda = cand.mean$lambda[max.idx], max = max.val))
  sd1.sel <- t(data.frame(alpha = cand.mean$alpha[sd1.idx], lambda = cand.mean$lambda[sd1.idx], sd1 = sd1.val))
  colnames(max.sel) <- c("ACC", "TPR(Hit)", "TPR(Miss)", "F score", "MCC", "AUC")
  colnames(sd1.sel) <- c("ACC", "TPR(Hit)", "TPR(Miss)", "F score", "MCC", "AUC")
  elapsed <- difftime(Sys.time(), start, units = "secs")[[1]]
  
  # Report
  if (verbose) {
    cat("\n\tHyperparameter Tuning Procedure for elastic net\n\n")
    cat(sprintf("  @data    : %s\n", deparse(substitute(data.mat))))
    cat(sprintf("  @fold    : %d\n", fold))
    cat(sprintf("  @iter    : %d\n", iter))
    cat(sprintf("  @parallel: %s\n", as.character(parallel)))
    cat(sprintf("  @seed    : %s\n", as.character(seed)))
    cat(sprintf("  @elapsed : %.4f (sec)\n", elapsed))
    
    cat("\n\tTuning Results (Maximum)\n\n")
    print(round(max.sel, 4))
    
    cat("\n\tTuning Results (SD1 rule)\n\n")
    print(round(sd1.sel, 4))
  }
  
  return(invisible(list(max.select = max.sel, sd1.select = sd1.sel, raw.mean = cand.mean, raw.sd = cand.sd, nzero = cand.nzero)))
}