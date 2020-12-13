\name{mfma}
\alias{mfma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting Mixtures of Factor Mixture Analyzers by the EM algorithm}
\description{
 \code{mfma} is used to estimate Mixtures of Factor Mixture Analyzers by the Expectation Maximization algorithm.
}
\usage{
mfma(y, k1, k2, r, it = 15, eps = 0.001, init = NULL,
    scaling = FALSE, seed = 7, model.names = c(1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A data matrix with \code{n} rows representing observations and \code{p} columns representing variables.}
  \item{k1}{The number of factor mixture analysis models}
  \item{k2}{The number of components of the factors of the factor mixture analysis model.}
  \item{r}{ The number of factors}
  \item{it}{The maximum number of iterations of the EM algorithm. By default it is set to 15.}
  \item{eps}{The lower bound for relative variation of the likelihood. It is used as alternative stopping rule for the EM algorithm: if the relative increment of the likelihood is lower than \code{eps} the fitting is stopped. The default is 0.0001.}
  \item{init}{A list containing initial values for all (of some) model parameters. If \code{NULL} (default) the algorithm starts from random values.}
  \item{scaling}{If TRUE (FALSE is default) the data are scaled before fitting the MFMA model.}
  \item{seed}{Fix the seed of the running. Default is 7.}
  \item{model.names}{ A vector describing the type of parameterization of the covariance matrices of the MFA and FMA parts. Default is c(1,1) which corresponds to
  the unrestricted covariances.  Other possibilities are: "2" is "EEE" and "EII" for the MFA and FMA covariances and "3" is VII.
  }
}

\value{A list containing the components:
  \item{h}{The number of estimated paramters}
  \item{k1}{The number of factor mixture analysis models}
  \item{k2}{The number of components of the factors of the factor mixture analysis model.}
  \item{H}{The estimated factor loading matrices}
  \item{lik}{The log-likelihood computed at each iteration of the EM algorithm}
  \item{w1}{The estimated weights of the MFA.}
  \item{w2}{The estimated weights of the FMA.}
  \item{muf}{A matrix with the estimated component means of MFA.}
  \item{mu}{A matrix with the estimated component means of FMA.}
  \item{sigma}{An array which contains the estimated component covariance of FMA.}
  \item{psi}{The noise diagonal variance matrix of MFA}
  \item{p}{The number of observed variables}
  \item{numobs}{The number of observations}
  \item{ph.y}{The posterior distribution of each mixture components.}
  \item{index}{The allocation vector.}
  \item{r}{The number of factors}
  \item{bic}{The BIC value}
  \item{aic}{The AIC value}
  \item{clc}{The CLC value}
  \item{icl.bic}{The ICL BIC value.}
  \item{s1,s2,s12}{Classification of the units performed through the three posteriors related to the double mixture}
  \item{cptime}{The computational time.}
}



\references{ C. Viroli, Dimensionally reduced model-based clustering through Mixtures of Factor Mixture Analyzers,
Journal of Classification, 27, 363-388 (2010).}
\author{Cinzia Viroli}
\keyword{multivariate}
