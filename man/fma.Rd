\name{fma}
\alias{fma}
\title{Fitting Factor Mixture Analysis model by the EM algorithm}
\description{
 \code{fma} is used to perform Factor Mixture Analysis (with covariates) on a matrix of data by the Expectation Maximization algorithm.
}
\usage{
fma(y, k, r, x.z = NULL, x.w = NULL, it = 15,
    eps = 1e-04, seed = 4, scaling = FALSE, init = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{y}{A data matrix with \code{n} rows representing observations and \code{p} columns representing variables.}
 \item{k}{The number of the mixture components.}
  \item{r}{The number of factors.}
  \item{x.z}{A matrix of covariates with \code{n} rows representing observations and \code{q} columns representing the predictors. These covariates are assumed to linearly affect the factor means.}
  \item{x.w}{A matrix of covariates with \code{n} rows representing observations and \code{m} columns representing the predictors. These covariates are assumed to differently affect the a priori probability of group membership.}
  \item{it}{The maximum number of iterations of the EM algorithm. By default it is set to 15.}
  \item{eps}{The lower bound for relative variation of the likelihood. It is used as alternative stopping rule for the EM algorithm: if the relative increment of the likelihood is lower than \code{eps} the fitting is stopped. The default is 0.0001.}
  \item{seed}{Fix the seed of the running. Default is 4.}
  \item{scaling}{If TRUE (FALSE is default) the data are scaled before fitting the FMA model.}
  \item{init}{A list containing initial values for all (of some) model parameters. If \code{NULL} (default) the algorithm starts from random values.}
}
\details{
Factor Mixture Analysis is a particular factor model with non Gaussian factors modelled by a multivariate Gaussian mixture. The \code{p} observed variables \code{y} are modelled in terms of the smaller set of \code{r} factors, \code{z}, and an additive specific term \code{u}: \code{y=Hz+u}, where \code{u} is assumed to be normally distributed with diagonal variance matrix \code{Psi}. \code{H} is the factor loading matrix. The model is fitted by the EM algorithm. The code implements also factor mixture model with covariates. Covariates may affect the observed variables into two manners: they are assumed to linearly affect the factor means (\code{x.z}) and \ or they can differently affect the a priori probability of group membership (\code{x.w}). The default is NULL which means that covariates are not incorporated in the model.}
\value{
\item{H}{The estimated factor loading matrix.}
  \item{lik}{The log-likelihood computed at each iteration of the EM algorithm.}
  \item{w}{A matrix with the estimated weights of the mixture.}
  \item{Beta}{An array of dimension \code{k x r x (q +1)} containing the vectors of regression
coefficients which are allowed to vary across the components. When no covariate is incorporated into the model,
i.e. \code{x.z} is NULL, \code{Beta} has dimension \code{k x r} and it correponds to the estimated component means of the mixture.}
  \item{phi}{A matrix of dimension \code{k x m} which contains the coefficients of the covariates affecting the group membership.}
  \item{sigma}{An array of dimension \code{k x r x r} which contains the estimated component covariance of the mixture.}
  \item{psi}{The noise diagonal variance matrix.}
  \item{ph.y}{The posterior distribution of each mixture components.}
  \item{z}{The reconstructed factor scores by Thomposon method.}
  \item{s}{The allocation vector.}
  \item{bic}{The BIC value.}
  \item{aic}{The AIC value.}
  \item{pc.time}{Computational time in seconds.}
}

\references{A. Montanari and C. Viroli (2010), Heteroscedastic Factor Mixture Analysis, Statistical Modelling, 10(4), 441-460.


A. Montanari and C. Viroli (2011), Dimensionally reduced mixtures of regression models, Journal of Statistical Planning and Inference, 141, 1744-1752.

            }

\author{Cinzia Viroli}
\examples{
library(MASS)
library(mvtnorm)
data(crabs)
y=as.matrix(crabs[,4:8])
fit=fma(y,k=2,r=2,it=50,eps=0.0001,scaling=TRUE) # without covariates
misc(fit$s,crabs[,1]) # compute the misclassification error

sex=as.double(crabs[,2])-1
fit=fma(y,k=2,r=2,x.z=sex,x.w=sex,it=50,eps=0.0001,scaling=TRUE) # with covariates
}
\keyword{multivariate}
