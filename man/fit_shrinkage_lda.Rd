\name{fit_shrinkage_lda}
\alias{fit_shrinkage_lda}
\title{Fit a shrinkage linear discriminant model}
\usage{
fit_shrinkage_lda(x, y, shrinkage = 0.1, priors = NULL)
}
\arguments{
\item{x}{Numeric matrix with samples in rows and features in columns.}
\item{y}{Class label for each row of \code{x}.}
\item{shrinkage}{Shrinkage fraction toward a diagonal within-class covariance
matrix, between 0 and 1.}
\item{priors}{Optional named or class-ordered class prior probabilities.}
}
\value{
An object containing class means, covariance, coefficients, intercepts,
priors, and the applied shrinkage.
}
\description{
Fits a linear discriminant model using a pooled within-class covariance
matrix shrunk toward its diagonal target.
}
\details{
Predicted class scores include the class-specific intercept and log prior;
they are not probabilities.
}
