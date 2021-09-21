library(deamer)
library(parallel)
library(doSNOW)
library(progress)

est.pow <- function(rand.gen, test, m, n, diff, iter = 1000, alpha = 0.05, seed = T, parallel = T, verbose = T) {
  # Setting up
  if (verbose) message("Setting up for simulation")
  
  if (seed) set.seed(0)
  
  start <- Sys.time()
  pow <- data.frame(diff = diff)
  
  if (parallel) { # Do parallel
    cl <- makeCluster(detectCores())
    registerDoSNOW(cl)
    
    if (seed) {
      seed.list <- round(1000 * rnorm(length(diff)))
    }
    
    prog <- function(n) if (verbose) setTxtProgressBar(txtProgressBar(max = length(diff), style = 3), n)
    
    if (verbose) message("Simulating")
    
    pow$pow <- foreach (i = seq_along(diff), .combine = c, .options.snow = list(progress = prog), .export = "depend.rand") %dopar% {
      if (seed) set.seed(seed.list[i])
      
      # Sampling
      samp <- matrix(rand.gen((m + n) * iter), nrow = iter, byrow = T)
      samp[,(m + 1):(m + n)] <- samp[,(m + 1):(m + n)] + diff[i]
      
      # Test
      test.res <- apply(samp, MARGIN = 1, function(x) test(x[1:m], x[(m + 1):(m + n)])$p.value)
      
      return(sum(test.res < alpha) / iter)
    }
    
    if (verbose) cat("\n")
    stopCluster(cl)
  } else {  # Do sequential
    pow$pow <- vector(length = length(diff))
    
    if (verbose) {
      message("Simulating")
      pb <- progress_bar$new(total = length(diff), clear = F)
    }
    
    for (i in seq_along(diff)) {
      if (verbose) pb$tick()
      
      # Sampling
      samp <- matrix(rand.gen((m + n) * iter), nrow = iter, byrow = T)
      samp[,(m + 1):(m + n)] <- samp[,(m + 1):(m + n)] + diff[i]
      
      # Test
      test.res <- apply(samp, MARGIN = 1, function(x) test(x[1:m], x[(m + 1):(m + n)])$p.value)
      pow$pow[i] <- sum(test.res < alpha) / iter
    }
  }
  
  elapsed <- difftime(Sys.time(), start, units = "secs")[[1]]
  
  # Report
  if (verbose) {
    cat("\n\tPower Estimation Procedure\n\n")
    cat(sprintf("  @iter   : %d\n", iter))
    cat(sprintf("  @m      : %d\n", m))
    cat(sprintf("  @n      : %d\n", n))
    cat(sprintf("  @seed   : %s\n", as.character(seed)))
    cat(sprintf("  @elapsed: %.4f (sec)\n", elapsed))
  }
  
  return(invisible(pow))
}

depend.rand <- function(n, g.size) {
  return(rep(rnorm(n / g.size), each = g.size) + 0.2 * rnorm(n))
}
