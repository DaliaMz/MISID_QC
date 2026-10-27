#' Fit a shrinkage linear discriminant model
#'
#' @param x Numeric matrix with samples in rows and features in columns.
#' @param y Class label for each row of `x`.
#' @param shrinkage Shrinkage fraction toward a diagonal within-class
#'   covariance matrix, between 0 and 1.
#' @param priors Optional named or class-ordered class prior probabilities.
#' @return An object containing class means, covariance, coefficients,
#'   intercepts, priors, and the applied shrinkage.
#' @export
fit_shrinkage_lda <- function(x, y, shrinkage = 0.1, priors = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  y <- as.character(y)

  if (length(dim(x)) != 2L || nrow(x) != length(y) || ncol(x) < 1L) {
    stop("x must be a non-empty feature matrix with one row per y value.")
  }
  if (any(!is.finite(x)) || anyNA(y) || any(!nzchar(y))) {
    stop("x and y must contain only finite, non-missing values.")
  }
  if (length(shrinkage) != 1L || !is.finite(shrinkage) ||
      shrinkage < 0 || shrinkage > 1) {
    stop("shrinkage must be a single number between 0 and 1.")
  }

  classes <- sort(unique(y))
  if (length(classes) < 2L) {
    stop("At least two classes are required.")
  }
  counts <- table(factor(y, levels = classes))
  if (nrow(x) <= length(classes)) {
    stop("More training samples than classes are required to estimate covariance.")
  }

  means <- do.call(rbind, lapply(classes, function(class) {
    colMeans(x[y == class, , drop = FALSE])
  }))
  rownames(means) <- classes
  residuals <- x - means[match(y, classes), , drop = FALSE]
  covariance <- crossprod(residuals) / (nrow(x) - length(classes))
  target <- diag(diag(covariance), nrow = ncol(x))
  covariance <- (1 - shrinkage) * covariance + shrinkage * target

  scale <- max(mean(diag(covariance)), 1)
  diag(covariance) <- diag(covariance) + scale * 1e-8
  inverse <- tryCatch(
    solve(covariance),
    error = function(e) stop("Could not regularize the within-class covariance.")
  )

  if (is.null(priors)) {
    priors <- as.numeric(counts) / sum(counts)
  } else {
    if (!is.null(names(priors))) {
      priors <- priors[classes]
    }
    if (length(priors) != length(classes) || anyNA(priors) ||
        any(!is.finite(priors)) || any(priors <= 0)) {
      stop("priors must provide a positive finite value for every class.")
    }
    priors <- priors / sum(priors)
  }
  names(priors) <- classes

  coefficients <- t(inverse %*% t(means))
  intercepts <- -0.5 * rowSums((means %*% inverse) * means) + log(priors)
  names(intercepts) <- classes

  structure(
    list(
      classes = classes,
      means = means,
      covariance = covariance,
      beta = coefficients,
      intercepts = intercepts,
      priors = priors,
      shrinkage = shrinkage,
      feature_names = colnames(x)
    ),
    class = "shrinkage_lda"
  )
}

#' Score samples with a shrinkage linear discriminant model
#'
#' Scores include the class-specific intercept and log prior, in contrast to
#' linear marginal contributions alone.
#'
#' @param object A model returned by [fit_shrinkage_lda()].
#' @param newdata Numeric matrix with samples in rows and matching features in
#'   columns.
#' @param ... Unused.
#' @return A matrix of complete class discriminant scores.
#' @export
predict_shrinkage_lda <- function(object, newdata, ...) {
  if (!inherits(object, "shrinkage_lda")) {
    stop("object must be a shrinkage_lda model.")
  }
  x <- as.matrix(newdata)
  storage.mode(x) <- "double"
  if (ncol(x) != ncol(object$beta) || any(!is.finite(x))) {
    stop("newdata must have finite values and the same number of features.")
  }
  if (!is.null(object$feature_names) && !is.null(colnames(x))) {
    if (!setequal(colnames(x), object$feature_names)) {
      stop("newdata feature names do not match the fitted model.")
    }
    x <- x[, object$feature_names, drop = FALSE]
  }
  scores <- x %*% t(object$beta)
  sweep(scores, 2L, object$intercepts, FUN = "+")
}
