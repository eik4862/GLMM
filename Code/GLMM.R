library(lme4)
library(progress)
library(parallel)
library(doSNOW)
library(reshape)

source("./Code/PredBase.R")

glmm.data.mat <- function(src) {
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
  data <- data[,-3]
  
  return(data)
}

glmm.cv <- function(formula, data.mat, fold = 10, iter = 100, prob.cut = seq(0, 1, by = 0.01), x.out = seq(0, 1, by = 0.01), seed = T, parallel = T, verbose = T) {
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
                        .packages = "lme4", .export = c("acc", "tpr.hit", "tpr.miss", "fscore", "mcc", "auc"),
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
        fit <- glmer(formula, data = data.mat[-hold.idx,], family = binomial, 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
        
        # Test
        pred <- predict(fit, data.mat[hold.idx,], type = "response")
        confu.mat <- table(pred = factor(ifelse(pred < 0.5, "Hit", "Miss"), levels = c("Hit", "Miss")),
                           real = factor(data.mat[hold.idx,]$res, levels = c("Hit", "Miss")))
        
        # Compute stats
        stats.tmp$acc[j] <- acc(confu.mat)
        stats.tmp$tpr.hit[j] <- tpr.hit(confu.mat)
        stats.tmp$tpr.miss[j] <- tpr.miss(confu.mat)
        stats.tmp$fscore[j] <- fscore(confu.mat)
        stats.tmp$mcc[j] <- mcc(confu.mat)
        
        # Compute ROC
        for (k in seq_along(prob.cut)) {
          confu.mat <- table(pred = factor(ifelse(pred < prob.cut[k], "Hit", "Miss"), levels = c("Hit", "Miss")),
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
        fit <- glmer(formula, data = data.mat[-hold.idx,], family = binomial, 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
        
        # Test
        pred <- predict(fit, data.mat[hold.idx,], type = "response")
        confu.mat <- table(pred = factor(ifelse(pred < 0.5, "Hit", "Miss"), levels = c("Hit", "Miss")),
                           real = factor(data.mat[hold.idx,]$res, levels = c("Hit", "Miss")))
        
        # Compute stats
        stats.tmp$acc[j] <- acc(confu.mat)
        stats.tmp$tpr.hit[j] <- tpr.hit(confu.mat)
        stats.tmp$tpr.miss[j] <- tpr.miss(confu.mat)
        stats.tmp$fscore[j] <- fscore(confu.mat)
        stats.tmp$mcc[j] <- mcc(confu.mat)
        
        # Compute ROC
        for (k in seq_along(prob.cut)) {
          confu.mat <- table(pred = factor(ifelse(pred < prob.cut[k], "Hit", "Miss"), levels = c("Hit", "Miss")),
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
    cat("\n\tCross Validation Procedure for GLMM\n\n")
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

sign.code <- function(p.val) {
  return(ifelse(p.val < 0.001, "***", ifelse(p.val < 0.01, "**", ifelse(p.val < 0.05, "*", ifelse(p.val < 0.1, ".", "")))))
}

glmm.lrt <- function(data.mat, vars, parallel = T, verbose = T) {
  start <- Sys.time()
  test.res <- data.frame(vars = vars)
  
  if (parallel) { # Do parallel
    # Setting up
    if (verbose) message("Setting up for LRT")
    
    prog <- function(n) if (verbose) setTxtProgressBar(txtProgressBar(max = length(vars), style = 3), n)
    cl <- makeCluster(detectCores())
    registerDoSNOW(cl)
    
    # Fit full model
    if (verbose) message("Fitting full model")
    full <- glmer(paste0("res~", paste(vars, collapse = "+"), "+(1|subj/sess)"), data = data.mat, family = binomial)
    
    if (verbose) message("Performing LRT")
    
    iter.test.res <- foreach(i = seq_along(vars), .combine = rbind, .packages = "lme4", .export = "sign.code",
                        .options.snow = list(progtest.ress = prog)) %dopar% {
      # Fit reduce model
      reduce <- glmer(paste0("res~", paste(vars[-i], collapse = "+"), "+(1|subj/sess)"), data = data.mat, family = binomial)
      
      # ANOVA
      test.res <- anova(reduce, full)[2,c(6, 8)]
      
      return(c(test.res, sign.code(test.res[2])))
    }
    
    stopCluster(cl)
    test.res <- cbind(test.res, iter.test.res)
    colnames(test.res)[2:4] <- c("stat", "p.val", "sign.code")
    rownames(test.res) <- NULL
  } else {  # Do sequential
    if (verbose) {
      message("Setting up for LRT")
      pb <- progtest.ress_bar$new(total = length(vars), clear = F)
    }
    
    test.res$stat <- vector(length = length(vars))
    test.res$p.val <- vector(length = length(vars))
    test.res$sign.code <- vector(length = length(vars))
    
    # Fit full model
    if (verbose) message("Fitting full model")
    full <- glmer(paste0("res~", paste(vars, collapse = "+"), "+(1|subj/sess)"), data = data.mat, family = binomial,
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
    
    if (verbose) message("Performing LRT")
    
    for (i in seq_along(vars)) {
      if (verbose) pb$tick()
      
      # Fit reduce model
      reduce.spec <- paste0("res~", paste(vars[-i], collapse = "+"), "+(1|subj/sess)")
      reduce <- glmer(reduce.spec, data = data.mat, family = binomial,
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
      
      # ANOVA
      test.res[i, c(2, 3)] <- anova(reduce, full)[2,c(6, 8)]
      test.res$sign.code[i] <- sign.code(test.res$p.val[i])
    }
  }
  
  elapsed <- difftime(Sys.time(), start, units = "secs")[[1]]
  
  # Report
  if (verbose) {
    cat("\n\tLikelihood Ratio Test for GLMM Coefficients\n\n")
    cat(sprintf("  @data    : %s\n", deparse(substitute(data.mat))))
    cat(sprintf("  @variabls: %d\n", length(vars)))
    cat(sprintf("  @parallel: %s\n", as.character(parallel)))
    cat(sprintf("  @elapsed : %.4f (sec)\n", elapsed))
    
    cat("\n\tTest results\n\n")
    print(test.res)
  }
  
  return(invisible(test.res))
}