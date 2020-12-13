\name{mfa}
\alias{mfa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting Mixtures of Factor Analyzers by the EM algorithm}
\description{
 \code{mfa} is used to estimate Mixtures of Factor Analyzers by the Expectation Maximization algorithm.
}
\usage{
mfa(y, k, r, it = 15, eps = 0.001, init = NULL,
    scaling = FALSE, seed = 7, model.names = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A data matrix with \code{n} rows representing observations and \code{p} columns representing variables.}
  \item{k}{The number of groups}
  \item{r}{ The number of factors}
  \item{it}{The maximum number of iterations of the EM algorithm. By default it is set to 15.}
  \item{eps}{The lower bound for relative variation of the likelihood. It is used as alternative stopping rule for the EM algorithm: if the relative increment of the likelihood is lower than \code{eps} the fitting is stopped. The default is 0.0001.}
  \item{init}{A list containing initial values for all (of some) model parameters. If \code{NULL} (default) the algorithm starts from random values.}
  \item{scaling}{If TRUE (FALSE is default) the data are scaled before fitting the MFMA model.}
  \item{seed}{Fix the seed of the running. Default is 7.}
  \item{model.names}{ A vector describing the type of parameterization of the covariance matrices of the MFA. Default is 1 which corresponds to the unrestricted covariances.  Other possibilities are: "2" is "EEE" and "3" is VII.}
}

\value{A list containing the components:
  \item{h}{The number of estimated paramters}
  \item{k}{The number of groups}
  \item{H}{The estimated factor loading matrices}
  \item{lik}{The log-likelihood computed at each iteration of the EM algorithm}
  \item{w}{The estimated weights of the MFA}
  \item{mu}{A matrix with the estimated component means of MFA}
  \item{psi}{The noise diagonal variance matrix of MFA}
  \item{p}{The number of observed variables}
  \item{numobs}{The number of observations}
  \item{r}{The number of factors}
  \item{bic}{The BIC value}
  \item{aic}{The AIC value}
  \item{clc}{The CLC value}
  \item{icl.bic}{The ICL BIC value}
  \item{s}{Posterior classification}
  \item{ps.y}{Posterior probabilities}
  \item{cptime}{The computational time}
}


\references{McLachlan, G., D. Peel, and R. Bean (2003). Modelling high-dimensional
data by mixtures of factor analyzers. Computational Statistics & Data
Analysis 41 (3), 379 - 388.
}

\author{Cinzia Viroli}
\examples{
library(MASS)
library(mvtnorm)
data(crabs)
y=as.matrix(crabs[,4:8])
fit=mfa(y,k=2,r=1,it=50,eps=0.0001,scaling=TRUE)
misc(fit$s,crabs[,1]) # compute the misclassification error
}
\keyword{multivariate}
