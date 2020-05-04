# Functions in the file return sets of high dimensional variable selection
# fitted models ("path")
# all the functions return a list(est.b = est.b)
# est.b is a m * (p+1) matrix, where m is the number of fitted models. All have
# column names (include intercept)


## ------- Lasso and elastic net --------
lm.lasso.path.helper <- function(X, Y, intercept, alpha) {
  # the cv.glmnet needs to take at least 2 covariates
  if (ncol(X) < 2){
    return (list(est.b = rbind(lm.ols(X, Y, intercept)$est.b,
                               lm.null(X, Y, intercept)$est.b)))
  }

  # set colnames
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  # fit
  mod <- glmnet::glmnet(X, Y, intercept = intercept, alpha = alpha)
  est.b <- cbind(mod$a0, t(as.matrix(mod$beta)))
  est.b <- est.b[which(!is.na(rowSums(est.b))), , drop = FALSE]

  # update column names
  colnames(est.b)[1] <- "intercept"
  colnames(est.b)[-1] <- colnames(X)

  return(list(est.b = est.b))
}

lm.lasso.path <- function(X, Y, intercept) {
  return(lm.lasso.path.helper(X, Y, intercept, alpha = 1))
}

lm.elastic.half.path <- function(X, Y, intercept) {
  return(lm.lasso.path.helper(X, Y, intercept, alpha = 0.5))
}

## ------- scad and mcp --------

lm.ncvreg.path.helper <- function(X, Y, intercept, penalty.method) {
  # set colnames
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  X <- as.matrix(X)
  class(X) <- "matrix"

  # fit
  mod <- ncvreg::ncvreg(X, Y, intercept = intercept,
                        family = "gaussian", penalty = penalty.method)

  # betas
  est.b <- t(as.matrix(mod$beta))
  if (!intercept) {
    est.b[, 1] <- 0
  }
  est.b <- est.b[which(!is.na(rowSums(est.b))), , drop = FALSE]

  # update column names
  colnames(est.b)[1] <- "intercept"
  colnames(est.b)[-1] <- colnames(X)

  return(list(est.b = est.b))
}

lm.scad.path <- function(X, Y, intercept) {
  return(lm.ncvreg.path.helper(X, Y, intercept, "SCAD"))
}

lm.mcp.path <- function(X, Y, intercept) {
  return(lm.ncvreg.path.helper(X, Y, intercept, "MCP"))
}

## ------- relaxo --------

lm.relaxo.path <- function(X, Y, intercept) {
  if (ncol(X) < 3){
    return (list(est.b = rbind(lm.ols(X, Y, intercept)$est.b,
                               lm.null(X, Y, intercept)$est.b)))
  }

  # set colnames
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  # remove variables that are not moving
  p <- ncol(X)
  var.name <- colnames(X)
  X <- as.matrix(X)

  sigma <- sapply(1:p, function(i) stats::sd(X[, i], na.rm = TRUE))
  i <- which(sigma != 0)
  X.dropped <- X[, i, drop = FALSE]

  # rescale X and Y (required by the package relaxo)
  scaled.Y <- scale(Y)
  scaled.X <- scale(X.dropped)
  X.scale <- attr(scaled.X, "scaled:scale")

  relaxo.mod <- relaxo::relaxo(scaled.X, scaled.Y)
  k <- length(relaxo.mod$lambda)
  j <- which(relaxo.mod$lambda[1:(k - 1)] != relaxo.mod$lambda[2:k]) + 1
  j <- c(1, j)

  # betas
  est.b <- matrix(0, nrow = length(j), ncol = p) # no intercept at the moment
  est.b[, i] <- relaxo.mod$beta[j, ] * attr(scaled.Y, "scaled:scale") / X.scale
  colnames(est.b) <- var.name # not include b0 yet

  # intercept
  b0 <- rep(0, nrow(est.b))
  if (intercept) {
    b0 <- c()
    for (k in seq_len(nrow(est.b)))
      b0 <- c(b0, mean(Y - X %*% est.b[k, ]))
  }

  est.b <- cbind(b0, est.b)
  est.b <- est.b[which(!is.na(rowSums(est.b))), , drop = FALSE]

  # update column names
  colnames(est.b)[1] <- "intercept"

  return(list(est.b = est.b))
}
