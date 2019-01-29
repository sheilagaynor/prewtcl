KMeansSparseCluster.permute2 <- function (x, y, K = 2, nperms = 25,
                                          wbounds = NULL, silent = FALSE,
                                          nvals = 10,
                                          ws0=rep(1/sqrt(ncol(x)), ncol(x)))
{
  if (is.null(wbounds))
    wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x)) * 0.9),
                       len = nvals))
  permx <- list()
  nnonzerows <- NULL
  nws <- sum(ws0!=0)
  for (i in 1:nperms) {
    permx[[i]] <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) permx[[i]][, j] <- sample(x[, j])
  }
  tots <- NULL
  out <- KMeansSparseCluster2(x, K, wbounds = wbounds, silent = silent,
                              ws0=ws0)
  for (i in 1:length(out)) {
    nnonzerows <- c(nnonzerows, sum(out[[i]]$ws != 0))
    bcss <- sparcl:::GetWCSS(x, out[[i]]$Cs)$bcss.perfeature
    tots <- c(tots, sum(out[[i]]$ws * bcss))
  }
  permtots <- matrix(NA, nrow = length(wbounds), ncol = nperms)
  for (k in 1:nperms) {
    if (!silent)
      cat("Permutation ", k, "of ", nperms, fill = TRUE)
    perm.out <- KMeansSparseCluster2(permx[[k]], K, wbounds = wbounds,
                                     silent = silent, ws0=ws0)
    for (i in 1:length(perm.out)) {
      perm.bcss <- sparcl:::GetWCSS(permx[[k]], perm.out[[i]]$Cs)$bcss.perfeature
      permtots[i, k] <- sum(perm.out[[i]]$ws * perm.bcss)
    }
  }
  gaps <- (log(tots) - apply(log(permtots), 1, mean))
  out <- list(tots = tots, permtots = permtots, nnonzerows = nnonzerows,
              gaps = gaps, sdgaps = apply(log(permtots), 1, sd), wbounds = wbounds,
              bestw = wbounds[which.max(gaps)])
  if (!silent)
    cat(fill = TRUE)
  class(out) <- "kmeanssparseperm"
  return(out)
}