\name{sig.features}
\alias{sig.features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate which features are significant for weighting.
}
\description{
The preweighted method requires a set of weights, which can be calculated by determining the association of initial clusters with matrix features.
}
\usage{
sig.features(x, Cs, perm = FALSE, B = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An nxp data matrix. There are n observations and p features.
}
  \item{Cs}{
Cluster assignments for n observations.
}
  \item{perm}{
The number of permutations.
}
  \item{B}{
The number of samples taken from the matrix.
}
}
\references{
Gaynor, Sheila and Bair, Eric, "Identification of biologically relevant subtypes via preweighted sparse clustering" (2013). The University of North Carolina at Chapel Hill Department of Biostatistics Technical Report Series. Working Paper 32.
}
\author{
 Sheila Gaynor and Eric Bair
}
\examples{
# Generate data
set.seed(13)
x <- matrix(rnorm(50*300),ncol=300)
x[1:25,1:50] <- x[1:25,1:50]+1
x[1:10,51:100] <- x[1:10,51:100]+2
x[40:50,51:100] <- x[40:50,51:100]+2
x <- scale(x, TRUE, TRUE)

# Determine weights
sparse.perm <- KMeansSparseCluster.permute(x, K = 2, nperm=5)
sparse <- KMeansSparseCluster(x, K = 2, wbounds=sparse.perm$bestw)
signif <- sig.features(x,sparse[[1]]$Cs)
threshold = 1*(signif > (0.05/300))
weights = threshold / sqrt(sum(threshold))

# Determine tuning parameter
perm <- KMeansSparseCluster.permute2(x, yval, K=2, ws0=weights, nperms=5)
index <- which(perm$gaps == max(perm$gaps))
tuning <- perm$wbound[index]

# Cluster the data
sspcl.out <- KMeansSparseCluster2(x, K=2, ws0=weights, wbounds=tuning)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
