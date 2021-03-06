\name{clustEff}
\alias{clustEff}
\title{
Cluster Effects Algorithm
}
\description{
This function implements the algorithm to cluster curves of effects obtained from a
quantile regression (qrcm; Frumento and Bottai, 2015) in which
the coefficients are described by flexible parametric functions of
the order of the quantile. This algorithm can be also used for clustering of curves observed in time, as in
functional data analysis.
}
\usage{
clustEff(Beta, p, alpha, k, ask=FALSE, k.min=1, k.max=min(10, (ncol(Beta)-1)),
        cluster.effects=TRUE, Beta.lower=NULL, Beta.upper=NULL,
        step=c("both", "shape", "distance"), plot=TRUE,
        method=c("ward.D", "ward.D2", "single", "complete", "average",
                 "mcquitty", "median", "centroid"))
}
\arguments{
  \item{Beta}{
    A matrix \kbd{n} x \kbd{q}. \kbd{q} represents the number of curves to cluster and \kbd{n} is either the length of percentiles used in the quantile regression or the length of the time vector.
  }
  \item{p}{
    The percentiles used in the quantile regression or the vector of time.
  }
  \item{alpha}{
    It is the alpha-percentile used for computing the dissimilarity matrix. If not fixed, the algorithm choose alpha=.25 (cluster.effects=TRUE) or alpha=.5 (cluster.effects=FALSE).
  }
  \item{k}{
    If fixed, it represents the number of clusters.
  }
  \item{ask}{
    If TRUE, after plotting the dendrogram, the user make is own choice about how many cluster to use.
  }
  \item{k.min}{
    The minimum number of clusters to let the algorithm to choose the best.
  }
  \item{k.max}{
    The maximum number of clusters to let the algorithm to choose the best.
  }
  \item{cluster.effects}{
    If TRUE, it selects the framework (quantile regression or curves clustering) in which to apply the clustering algorithm. }
  \item{Beta.lower}{
    A matrix \kbd{n} x \kbd{q}. \kbd{q} represents the number of lower interval of the curves to cluster and \kbd{n} the length of percentiles used in quantile regression. Used only if cluster.effects=TRUE.
  }
  \item{Beta.upper}{
    A matrix \kbd{n} x \kbd{q}. \kbd{q} represents the number of upper interval of the curves to cluster and \kbd{n} the length of percentiles used in quantile regression. Used only if cluster.effects=TRUE.
  }
  \item{step}{
    The steps used in computing the dissimilarity matrix. Default is "both"=("shape" and  "distance")
  }
  \item{plot}{
    If TRUE, dendrogram, boxplot and clusters are plotted.
  }
  \item{method}{
    The agglomeration method to be used.
  }
}
\details{
  Quantile regression models conditional quantiles of a response variabile,
  given a set of covariates. Assume that each coefficient can be expressed as a parametric function of      \eqn{p} in the form: \deqn{\beta(p | \theta) = \theta_{0} + \theta_1 b_1(p) + \theta_2 b_2(p) +      \ldots}{ \beta(p | \theta) = \theta0 + \theta1*b1(p) + \theta2*b2(p) + \ldots}
  where \eqn{b_1(p), b_2(p, \ldots)}{b1(p), b2(p), \ldots} are known functions of \eqn{p}.
}
\value{
An object of class \dQuote{\code{clustEff}}, a list containing the following items:
\item{call}{the matched call.}
\item{X}{The curves matrix.}
\item{X.mean}{The mean curves matrix of dimension \kbd{n} x \kbd{k}.}
\item{X.mean.dist}{The within cluster distance from the mean curve.}
\item{X.lower}{The lower interval matrix.}
\item{X.mean.lower}{The mean lower interval of dimension \kbd{n} x \kbd{k}.}
\item{X.upper}{The upper interval matrix.}
\item{X.mean.upper}{The mean upper interval of dimension \kbd{n} x \kbd{k}.}
\item{k}{The number of selected clusters.}
\item{p}{The percentiles used in quantile regression coefficient modeling or the time otherwise.}
\item{diss.matrix}{The dissimilarity matrix.}
\item{oggSilhouette}{An object of class \dQuote{\code{silhouette}}.}
\item{oggHclust}{An object of class \dQuote{\code{hclust}}.}
\item{clusters}{The vector of clusters.}
\item{alpha.dist}{The vector of alpha-distances corresponding to the alpha-percentile of the distances along the percentiles.}
\item{distance}{A vector of goodness measures used to select the best number of clusters.}
\item{step}{The selected step.}
\item{method}{The used agglomeration method.}
\item{cut.method}{The used method to select the best number of clusters.}
\item{alpha}{The selected alpha-percentile.}
}
\references{
Sottile, G and Adelfio, G (2017). \emph{Clustering of effects through quantile regression}.
Proceedings: International Workshop of Statistical Modeling.

Frumento, P., and Bottai, M. (2015). \emph{Parametric modeling of quantile regression coefficient functions}.
Biometrics, doi: 10.1111/biom.12410.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{summary.clustEff}}, \code{\link{plot.clustEff}},
for summary and plotting.
\code{\link{extract.object}} to extract useful objects for the clustering algorithm through a quantile regression coefficient modeling in a multivariate case.
}
\examples{

  ##### Using simulated data in all examples

  # see the documentation for 'clustEff-package'

}
\keyword{models}
\keyword{clustering algorithm}
