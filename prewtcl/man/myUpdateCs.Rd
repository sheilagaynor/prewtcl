\name{myUpdateCs}
\alias{myUpdateCs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to update the cluster assignments.
}
\description{
Cluster assignments from kmeans updated.
}
\usage{
myUpdateCs(x, K, ws, nstart = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 An nxp data matrix. There are n observations and p features.
}
  \item{K}{
The number of clusters desired ("K" in K-means clustering).
}
  \item{ws}{
%%     ~~Describe \code{ws} here~~
}
  \item{nstart}{
%%     ~~Describe \code{nstart} here~~
}
}
\references{
Gaynor, Sheila and Bair, Eric, "Identification of biologically relevant subtypes via preweighted sparse clustering" (2013). The University of North Carolina at Chapel Hill Department of Biostatistics Technical Report Series. Working Paper 32.
}
\author{
 Sheila Gaynor and Eric Bair
}
\examples{
## The function is currently defined as
function (x, K, ws, nstart=10) 
{
    Cs <- kmeans(x, centers=K, nstart=nstart)$cluster
    UpdateCs(x, K, ws, Cs)
  }
  #UpdateCs is a hidden function of sparcl.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
