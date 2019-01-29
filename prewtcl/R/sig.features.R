sig.features <- function(x, Cs, perm=FALSE, B=1000)
{ library(sparcl)
  
  if (!perm) {
    n.Cs <- length(unique(Cs))
    junk <- sparcl:::GetWCSS(x, Cs)
    pvals <- 1-pf((junk$bcss.perfeature/(n.Cs-1))/
                    (junk$wcss.perfeature/(nrow(x)-n.Cs)), n.Cs-1,
                  nrow(x)-n.Cs)}
  else {
    x <- scale(x)
    x.bcss <- sparcl:::GetWCSS(x, Cs)$bcss.perfeature
    
    if (TRUE) {
      nf <- ncol(x)
      xstar.bcss <- matrix(nrow=B, ncol=nf)
      for (i in 1:B) {
        x.star <- matrix(nrow=nrow(x), ncol=ncol(x))
        for (j in 1:ncol(x)) {
          x.star[,j] <- sample(x[,j])
        }
        xstar.bcss[i,] <- sparcl:::GetWCSS(x.star, Cs)$bcss.perfeature
      }
      pvals <- rep(NA, nf)
      for (i in 1:nf) {
        pvals[i] <- sum(x.bcss[i]<=xstar.bcss)/(B*nf)
      }    }  }
  return(pvals)
}