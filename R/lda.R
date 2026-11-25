#' Fit a shrinkage linear discriminant model
#'
#' @param x Numeric matrix with samples in rows and features in columns.
#' @param y Class label for each row of `x`.
#' @param shrinkage Shrinkage fraction toward a diagonal within-class
#'   covariance matrix, between 0 and 1.
#' @param priors Optional named or class-ordered class prior probabilities.
#' @return An object containing class means, coefficients, intercepts, priors,
#'   and the regularized covariance diagonal. A dual-space calculation avoids
#'   constructing a features-by-features covariance matrix.
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
  degrees_of_freedom <- nrow(x) - length(classes)
  variance <- colSums(residuals^2) / degrees_of_freedom
  regularizer <- max(mean(variance), 1) * 1e-8
  covariance_diagonal <- shrinkage * variance + regularizer
  coefficients_transposed <- sweep(
    t(means), 1L, covariance_diagonal, FUN = "/"
  )
  middle <- NULL
  if (shrinkage < 1) {
    weighted_residuals <- sweep(
      residuals, 2L, covariance_diagonal, FUN = "/"
    )
    gram <- residuals %*% t(weighted_residuals)
    middle <- tryCatch(
      solve(diag(degrees_of_freedom / (1 - shrinkage),
                 nrow = nrow(x)) + gram),
      error = function(e) stop("Could not regularize the within-class covariance.")
    )
    coefficients_transposed <- coefficients_transposed -
      sweep(t(residuals), 1L, covariance_diagonal, FUN = "/") %*%
      middle %*%
      (residuals %*% coefficients_transposed)
  }

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

  coefficients <- t(coefficients_transposed)
  intercepts <- -0.5 * rowSums(coefficients * means) + log(priors)
  names(intercepts) <- classes

  structure(
    list(
      classes = classes,
      means = means,
      covariance_diagonal = covariance_diagonal,
      residuals = residuals,
      middle = middle,
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

# Distance from each row of `newdata` to the nearest class mean under the
# model's shrinkage covariance (Mahalanobis, on the original scale). With
# S = covariance, the squared distance to class k is
# x' S^-1 x - 2 (score_k - log prior_k), and x' S^-1 x follows from the same
# Woodbury identity used when fitting, so S is never formed.
.nearest_class_distance <- function(model, newdata) {
  x <- as.matrix(newdata)
  if (!is.null(model$feature_names) && !is.null(colnames(x))) {
    x <- x[, model$feature_names, drop = FALSE]
  }
  scores <- predict_shrinkage_lda(model, x)
  scaled <- sweep(x, 2L, model$covariance_diagonal, FUN = "/")
  quadratic <- rowSums(scaled * x)
  if (!is.null(model$middle)) {
    projected <- model$residuals %*% t(scaled)
    quadratic <- quadratic - colSums(projected * (model$middle %*% projected))
  }
  best <- apply(sweep(scores, 2L, log(model$priors), FUN = "-"), 1L, max)
  sqrt(pmax(quadratic - 2 * best, 0))
}

# Pooled Mahalanobis distance to the cohort center. Unlike
# `.nearest_class_distance()`, this statistic does not reward fitting an
# unknown sample to the nearest enrolled class.
.global_mahalanobis_distance <- function(model, newdata) {
  x <- as.matrix(newdata)
  if (!is.null(model$feature_names) && !is.null(colnames(x))) {
    x <- x[, model$feature_names, drop = FALSE]
  }
  center <- colSums(sweep(model$means, 1L, model$priors, FUN = "*"))
  residual <- sweep(x, 2L, center, FUN = "-")
  scaled <- sweep(residual, 2L, model$covariance_diagonal, FUN = "/")
  quadratic <- rowSums(scaled * residual)
  if (!is.null(model$middle)) {
    projected <- model$residuals %*% t(scaled)
    quadratic <- quadratic - colSums(projected * (model$middle %*% projected))
  }
  sqrt(pmax(quadratic, 0))
}

# Reconstruction residual after projecting pooled within-class variation onto
# its leading principal components. Large residuals represent variation outside
# the enrolled covariance subspace.
.pca_reconstruction_residual <- function(model, newdata, n_components = 10L) {
  x <- as.matrix(newdata)
  if (!is.null(model$feature_names) && !is.null(colnames(x))) {
    x <- x[, model$feature_names, drop = FALSE]
  }
  center <- colSums(sweep(model$means, 1L, model$priors, FUN = "*"))
  scale <- sqrt(model$covariance_diagonal)
  training <- sweep(model$residuals, 2L, scale, FUN = "/")
  decomposition <- svd(training, nu = 0L,
                       nv = min(as.integer(n_components),
                                nrow(training), ncol(training)))
  k <- ncol(decomposition$v)
  if (!k) return(rep(0, nrow(x)))
  standardized <- sweep(sweep(x, 2L, center, FUN = "-"), 2L, scale,
                        FUN = "/")
  coefficients <- standardized %*% decomposition$v
  residual <- standardized - coefficients %*% t(decomposition$v)
  sqrt(rowMeans(residual * residual))
}

# Mean distance to the k nearest training residuals after nearest-class
# centering and diagonal shrinkage-covariance standardization.
.nearest_training_distance <- function(model, newdata, k = 5L) {
  x <- as.matrix(newdata)
  if (!is.null(model$feature_names) && !is.null(colnames(x))) {
    x <- x[, model$feature_names, drop = FALSE]
  }
  scale <- sqrt(model$covariance_diagonal)
  class_distance <- do.call(cbind, lapply(seq_len(nrow(model$means)), function(i) {
    residual <- sweep(x, 2L, model$means[i, ], FUN = "-")
    rowSums((residual / scale) * (residual / scale))
  }))
  nearest <- max.col(-class_distance, ties.method = "first")
  nearest_means <- t(vapply(seq_len(nrow(x)), function(i) {
    model$means[nearest[i], ]
  }, numeric(ncol(model$means))))
  test_residual <- x - nearest_means
  test_residual <- sweep(test_residual, 2L, scale, FUN = "/")
  training_residual <- sweep(model$residuals, 2L, scale, FUN = "/")
  squared <- rowSums(test_residual * test_residual)
  squared <- sweep(-2 * (test_residual %*% t(training_residual)),
                   1L, squared, FUN = "+")
  squared <- sweep(squared, 2L, rowSums(training_residual * training_residual),
                   FUN = "+")
  squared <- pmax(squared, 0)
  k <- min(as.integer(k), nrow(training_residual))
  sqrt(apply(squared, 1L, function(row) mean(sort(row)[seq_len(k)])))
}
