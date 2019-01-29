\name{KMeansSparseCluster2}
\alias{KMeansSparseCluster2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs preweighted sparse k-means clustering
}
\description{
This function performs sparse k-means clustering. The user must specify the number of clusters K, an L1 	bound on w (which are the feature weights), and the initial set of weights.
}
\usage{
KMeansSparseCluster2(x, K, wbounds = NULL, nstart = 20, silent = FALSE, maxiter = 6, ws0 = rep(1/sqrt(ncol(x)), ncol(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     An nxp data matrix. There are n observations and p features.
}
  \item{K}{
    The number of clusters desired ("K" in K-means clustering).
}
  \item{wbounds}{
    A single L1 bound on w (the feature weights), or a vector of L1 bounds on w.
	If wbound is small, then few features will have non-zero weights. If wbound is
	large then all features will have non-zero weights. Should be greater than 1.
}
  \item{nstart}{
     The number of random starts for the k-means algorithm.
}
  \item{silent}{
    Whether or not to print out progress.
}
  \item{maxiter}{
     The maximum number of iterations.
}
  \item{ws0}{
     The set of initial weights to be used in clustering.
}
}
\value{
  For each value of wbounds, the elements are:
  \item{ws }{The p-vector of the feature weights.}
  \item{comp2 }{The clustering assignments obtained.}

}
\references{
Gaynor, Sheila and Bair, Eric, "Identification of biologically relevant subtypes via preweighted sparse clustering" (2013). The University of North Carolina at Chapel Hill Department of Biostatistics Technical Report Series. Working Paper 32.
}
\author{
 Sheila Gaynor and Eric Bair
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{KMeansSparseCluster.permute2}}
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
