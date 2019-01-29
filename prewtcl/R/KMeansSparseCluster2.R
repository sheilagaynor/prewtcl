KMeansSparseCluster2 <-
  function (x, K, wbounds = NULL, nstart = 20, silent = FALSE, 
            maxiter = 6, ws0 = rep(1/sqrt(ncol(x)), ncol(x))) 
  {
    if (is.null(wbounds)) 
      wbounds <- seq(1.5, sqrt(ncol(x)), len = 20)
    if (min(wbounds) < 1) 
      stop("wbounds should be greater than or equal to 1")
    wbounds <- c(wbounds)
    out <- list()
    Cs <- kmeans(x, centers = K, nstart = nstart)$cluster
    for (i in 1:length(wbounds)) {
      if (length(wbounds) > 1 && !silent) 
        cat(i, fill = FALSE)
      ws <- ws0
      ws.old <- rnorm(ncol(x))
      store.bcss.ws <- NULL
      niter <- 0
      while ((sum(abs(ws - ws.old))/sum(abs(ws.old))) > 1e-04 && 
               niter < maxiter) {
        if (!silent) 
          cat(niter, fill = FALSE)
        niter <- niter + 1
        ws.old <- ws
        Cs <- sparcl:::UpdateCs(x, K, ws, Cs)
        ws <- sparcl:::UpdateWs(x, Cs, wbounds[i])
        store.bcss.ws <- c(store.bcss.ws, sum(sparcl:::GetWCSS(x, 
                                                      Cs)$bcss.perfeature * ws))
      }
      out[[i]] <- list(ws = ws, Cs = Cs, wcss = sparcl:::GetWCSS(x, 
                                                        Cs, ws), crit = store.bcss.ws, wbound = wbounds[i])
    }
    if (!silent) 
      cat(fill = TRUE)
    if (length(wbounds) == 1) {
      out <- out[[1]]
      class(out) <- "kmeanssparse"
      return(out)
    }
    class(out) <- "multikmeanssparse"
    return(out)
  }
