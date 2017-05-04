\name{clustEff-package}
\alias{clustEff-package}
\docType{package}
\title{
Clusters of effects
}
\description{
This package implements a general algorithm to cluster curves of effects obtained from a
quantile regression (qrcm; Frumento and Bottai, 2015) in which
the coefficients are described by flexible parametric functions of
the order of the quantile. This algorithm can be also used for clustering of curves observed in time, as in
functional data analysis.
}
\details{
\tabular{ll}{
Package: \tab clustEff\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-03-06\cr
License: \tab GPL-2\cr
}
The function \code{\link{clustEff}} allows to specify the type of the curves to apply the proposed clustering algorithm. The function \code{\link{extract.object}}  extracts the matrices, in case of multivariate response, through the quantile regression coefficient modeling, useful to run the main algorithm.  The auxiliary functions \code{\link{summary.clustEff}} and \code{\link{plot.clustEff}}
can be used to extract information from the main algorithm.
}
\author{
Gianluca Sottile

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
\references{
Sottile, G and Adelfio, G (2017). \emph{Clustering of effects through quantile regression}.
Proceedings International Workshop of Statistical Modeling.

Frumento, P., and Bottai, M. (2015). \emph{Parametric modeling of quantile regression coefficient functions}.
Biometrics, doi: 10.1111/biom.12410.}
\keyword{ package }
\examples{

# use simulated data
# CURVES EFFECTS CLUSTERING

set.seed(1234)
n <- 1000
q <- 2
k <- 5
x1 <- runif(n, 0, 5)
x2 <- runif(n, 0, 5)

X <- cbind(x1, x2)
rownames(X) <- 1:n
colnames(X) <- paste0("X", 1:q)

theta1 <- matrix(c(1, 1, 0, 0, 0, .5, 0, .5, 1, 2, .5, 0, 2, 1, .5),
                 ncol=k, byrow=TRUE)

theta2 <- matrix(c(1, 1, 0, 0, 0, -.3, 0, .5, 1, .5, -1.5, 0, -1, -.5, 1),
                 ncol=k, byrow=TRUE)

theta3 <- matrix(c(1, 1, 0, 0, 0, .3, 0, -.5, -1, 2, -.5, 0, 1, -.5, -1),
                 ncol=k, byrow=TRUE)

rownames(theta3) <- rownames(theta2) <- rownames(theta1) <-
    c("(intercept)", paste("X", 1:q, sep=""))
colnames(theta3) <- colnames(theta2) <- colnames(theta1) <-
    c("(intercept)", "qnorm(p)", "p", "p^2", "p^3")

Theta <- list(theta1, theta2, theta3)

B <- function(p, k){matrix(cbind(1, qnorm(p), p, p^2, p^3), nrow=k, byrow=TRUE)}
Q <- function(p, theta, B, k, X){rowSums(X * t(theta \%*\% B(p, k)))}

s <- matrix(1, q+1, k)
s[2:(q+1), 2] <- 0
s[1, 3:k] <- 0

Y <- matrix(NA, nrow(X), 15)
for(i in 1:15){
  if(i <= 5) Y[, i] <- Q(runif(n), Theta[[1]], B, k, cbind(1, X))
  if(i <= 10 & i > 5) Y[, i] <- Q(runif(n), Theta[[2]], B, k, cbind(1, X))
  if(i <= 15 & i > 10) Y[, i] <- Q(runif(n), Theta[[3]], B, k, cbind(1, X))
}

XX <- extract.object(Y, X, intercept=TRUE, formula.p= ~ I(p) + I(p^2) + I(p^3))
seqP <- XX$p

obj <- clustEff(XX$X$X1, seqP, Beta.lower=XX$Xl$X1, Beta.upper=XX$Xr$X1)
summary(obj)
plot(obj, xvar="clusters", add=TRUE)
plot(obj, xvar="clusters", add=FALSE)
plot(obj, xvar="dendrogram")
plot(obj, xvar="boxplot")

# CURVES CLUSTERING IN FUNCTIONAL DATA ANALYSIS

set.seed(1234)
n <- 1000
x <- 1:n/n

Y <- matrix(0, n, 30)

sigma2 <- 4*pmax(x-.2, 0) - 8*pmax(x-.5, 0) + 4*pmax(x-.8, 0)

mu <- sin(3*pi*x)
for(i in 1:10) Y[, i] <- mu + rnorm(length(x), 0, pmax(sigma2, 0))

mu <- cos(3*pi*x)
for(i in 11:23) Y[,i] <- mu + rnorm(length(x), 0, pmax(sigma2,0))

mu <- sin(3*pi*x)*cos(pi*x)
for(i in 24:28) Y[, i] <- mu + rnorm(length(x), 0, pmax(sigma2, 0))

mu <- 0 #sin(1/3*pi*x)*cos(2*pi*x)
for(i in 29:30) Y[, i] <- mu + rnorm(length(x), 0, pmax(sigma2, 0))

obj2 <- clustEff(Y, x, cluster.effects=FALSE)
summary(obj2)

plot(obj2, xvar="clusters", add=TRUE)
plot(obj2, xvar="clusters", add=FALSE)
plot(obj2, xvar="dendrogram")
plot(obj2, xvar="boxplot")

}



