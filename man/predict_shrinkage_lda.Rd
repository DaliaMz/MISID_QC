\name{predict_shrinkage_lda}
\alias{predict_shrinkage_lda}
\title{Score samples with a shrinkage linear discriminant model}
\usage{
predict_shrinkage_lda(object, newdata, ...)
}
\arguments{
\item{object}{A model returned by \code{\link{fit_shrinkage_lda}}.}
\item{newdata}{Numeric matrix with samples in rows and matching features in
columns.}
\item{...}{Unused.}
}
\value{
A matrix of complete class discriminant scores, including intercepts and
class priors.
}
\description{
Calculates a discriminant score for every input sample and fitted class.
}
