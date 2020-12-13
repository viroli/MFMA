\name{misc}
\alias{misc}
\title{Misclassification error}
\description{
  An internal function which computes the misclassification error between two partitions
}
\usage{
misc(classification, truth)
}
\arguments{
  \item{classification}{A numeric or character vector of class labels.}
  \item{truth}{A numeric or character vector of truth class labels. The length of truth should be the same as that of classification.}
}
\value{
The misclassification error (a scalar).
}
\seealso{\code{\link{fma}}}
\keyword{multivariate}
