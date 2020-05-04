#' Get fitted models by fitting some variable selection methods
#' @export lm.compare.method
#' @param X covariates (n times p matrix, n: number of entries, p: number of covariates)
#' @param Y response (vector with n entries)
#' @param intercept TRUE to fit the data with an intercept, FALSE to fit the data without an intercept
#' @param method.names vector of method names to be used for fitting. Choose among "lasso", "elastic", "relaxo", "mcp" and "scad". Default is to fit the data using all methods listed above
#' @param log.level log level to set. Default is NULL, which means no change in log level. See the function CSUV::set.log.level for more details
#' @return estimated coefficients in a form of matrix. Each row corresponds to a method and each column corresponds to a covariate, with the first column corresponds to the intercept
#' @examples
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' compare.mod = lm.compare.method(X, Y, intercept = FALSE)
#' print(compare.mod)
lm.compare.method <- function(X, Y, intercept, method.names = NULL, log.level = NULL) {
  if (!is.null(log.level)){
    set.log.level(log.level)
  }
  if (is.null(method.names)) {
    method.names <- names(get.compare.methods())
  }
  return(get.compare.fit(X, Y, intercept, method.names, current.compare.fit = NULL))
}

lm.null <- function(X, Y, intercept) {
  # return: est.b of length p+1
  ols.mod <- NULL
  est.b <- double()

  if (is.null(X) || !ncol(X)) {
    if (intercept) {
      ols.mod <- stats::lm(Y ~ 1) # intercept only
      est.b <- ols.mod$coefficients
    } else{
      est.b <- 0
    }
  } else{
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("X", seq_len(ncol(X)))
    }
    d <- data.frame(X, Y)
    if (intercept) {
      ols.mod <- stats::lm(Y ~ 1, data = d)
      est.b <- c(ols.mod$coefficients, rep(0, ncol(X)))
    } else{
      ols.mod <- stats::lm(Y ~ 0, data = d)
      est.b <- rep(0, ncol(X)+1)
    }
  }
  value <- list(raw.mod = ols.mod, est.b = est.b)
  attr(value, "class") <- "lm.result"
  return(value)
}

lm.ols <- function(X, Y, intercept) {
  # return: est.b of length p+1
  ols.mod <- NULL
  est.b <- double()

  if (is.null(X) || !ncol(X)) {
    if (intercept) {
      ols.mod <- stats::lm(Y ~ 1) # intercept only
      est.b <- ols.mod$coefficients
    } else{
      est.b <- 0
    }
  } else{
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("X", seq_len(ncol(X)))
    }
    d <- data.frame(X, Y)
    if (intercept) {
      ols.mod <- stats::lm(Y ~ ., data = d)
      est.b <- ols.mod$coefficients
    } else{
      ols.mod <- stats::lm(Y ~ . + 0, data = d)
      est.b <- c(0, ols.mod$coefficients)
    }
    est.b[which(is.na(est.b))] <- 0
  }
  value <- list(raw.mod = ols.mod, est.b = est.b)
  attr(value, "class") <- "lm.result"
  return(value)
}

## ======== variable selection methods =========

lm.lasso.min <- function(X, Y, intercept) {
  # for some unknown reason, the cv.glmnet needs to take at least 2 covariates
  if (ncol(X) < 2)
    return(lm.ols(X, Y, intercept))

  # set colnames
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  # fit
  lasso.mod <- glmnet::cv.glmnet(X, Y, intercept = intercept)

  # betas
  est.b <- as.numeric(glmnet::coef.glmnet(lasso.mod, s = "lambda.min"))
  names(est.b[-1]) <- colnames(X)

  value <- list(raw.mod = lasso.mod, est.b = est.b)
  attr(value, "class") <- "lm.result"
  return(value)
}


lm.lasso.adapt <- function(X, Y, intercept) {
  if (ncol(X) < 2)
    return(lm.ols(X, Y, intercept))

  # set colnames
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  adapt.mod <- parcor::adalasso(X = X, y = Y, intercept = intercept)
  est.b <- c(adapt.mod$intercept.adalasso,
             adapt.mod$coefficients.adalasso)
  names(est.b[-1]) <- colnames(X)

  value <- list(raw.mod = adapt.mod, est.b = est.b)
  attr(value, "class") <- "lm.result"
  return(value)
}

lm.elastic.half <- function(X, Y, intercept) {

  # for some unknown reason, the cv.glmnet needs to take at least 2 covariates
  if (ncol(X) < 2)
    return(lm.ols(X, Y, intercept))

  # set colnames
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

  # fit
  lasso.mod <- glmnet::cv.glmnet(X, Y, intercept = intercept, alpha = 0.5)

  # betas
  est.b <- as.numeric(glmnet::coef.glmnet(lasso.mod, s = "lambda.min"))
  names(est.b[-1]) <- colnames(X)

  value <- list(raw.mod = lasso.mod, est.b = est.b)
  attr(value, "class") <- "lm.result"
  return(value)
}

lm.ncvreg.helper <- function(X, Y, intercept, penalty.method) {
  # set colnames
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

  X <- as.matrix(X)
  class(X) <- "matrix"

  # fit
  mod <- ncvreg::cv.ncvreg(X, Y, intercept = intercept,
                           family = "gaussian", penalty = penalty.method)

  # betas
  est.b <- mod$fit$beta[, mod$min]
  if (!intercept) est.b[1] <- 0

  value <- list(raw.mod = mod, est.b = est.b)
  attr(value, "class") <- "lm.result"
  return(value)
}

lm.scad <- function(X, Y, intercept) {
  return(lm.ncvreg.helper(X, Y, intercept, "SCAD"))
}

lm.mcp <- function(X, Y, intercept) {
  return(lm.ncvreg.helper(X, Y, intercept, "MCP"))
}

lm.relaxo <- function(X, Y, intercept) {
  # relaxo does not work for < 3 covariates
  if (ncol(X) < 3)
    return(lm.ols(X, Y, intercept))

  # set colnames
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

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

  relaxo.mod <- cvrelaxo1(scaled.X, scaled.Y, warn = FALSE)

  # betas
  est.b <- rep(0, p) # no intercept at the moment
  est.b[i] <- as.numeric(relaxo.mod$beta * attr(scaled.Y, "scaled:scale") / X.scale)
  names(est.b) <- var.name # not include b0 yet

  # intercept
  b0 <- 0
  if (intercept) {
    b0 <- mean(Y - X %*% est.b)
  }

  value <- list(raw.mod = relaxo.mod, est.b = c(b0, est.b))
  attr(value, "class") <- "lm.result"
  return(value)
}

## ======== refit =========
#' Get the ordinary least square estimated coefficients on a set of previously selected covariates
#' @export lm.ols.refit
#' @param X covariates (n times p matrix, n: number of entries, p: number of covariates)
#' @param Y response (vector with n entries)
#' @param intercept TRUE to fit the data with an intercept, FALSE to fit the data without an intercept
#' @param est.betas estimated betas from previous fitted result. It can be a vector with p+1 entries (first entry as intercept) or a matrix with p+1 columns. Non-zero coefficient means the corresponding covariate is selected
#' @param log.level log level to set. Default is NULL, which means no change in log level. See the function CSUV::set.log.level for more details
#' @return a list of estimated coefficients
#' @examples
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' est.beta = rep(0, 11)
#' est.beta[2:5] = 1
#' ols.mod = lm.ols.refit(X, Y, intercept = FALSE, est.betas = est.beta)
#' print(ols.mod$est.b)
lm.ols.refit <- function(X, Y, intercept, est.betas, log.level = NULL) {
  if (!is.null(log.level)){
    set.log.level(log.level)
  }
  if (is.vector(est.betas)) {
    if (ncol(X) != length(est.betas) - 1) {
      stop("wrong number of variables for lm.ols.refit")
    }
    return(lm.ols.refit.one(X, Y, intercept, est.betas))
  } else{
    if (ncol(X) != ncol(est.betas) - 1) {
      stop("wrong number of variables for lm.ols.refit")
    }
    result <- do.call(rbind, lapply(seq_len(nrow(est.betas)), function(i)
      lm.ols.refit.one(X, Y, intercept, est.betas[i, ])$est.b))
    return(list(est.b = result))
  }
}

lm.ols.refit.one <- function(X, Y, intercept, est.b) {
  num.var.sel <- sum(est.b[-1] != 0)
  if (num.var.sel < length(Y)) { # because ols can only handle p<n
    if (length(est.b[-1]) != ncol(X)) {
      # print(est.b[-1])
      # print(ncol(X))
      stop("wrong number of variables for lm.ols.refit.one")
    }
    i <- which(est.b[-1] != 0)
    update.index <- 1 # intercept
    if (length(i)) {
      update.index <- c(update.index, (i + 1))
    }
    ols.mod <- lm.ols(X[, i, drop = FALSE], Y, intercept)
    est.b[update.index] <- ols.mod$est.b

    return(list(raw.mod = ols.mod$raw.mod, est.b = est.b))
  } else{
    warning(paste0("cannot refit as the number of observations (", length(Y), "is not greater than the number of covariates", num.var.sel, ". Use the est beta instead"))
    return(list(est.b = est.b))
  }
}

## ======== mse and prediction =========
#' Calculate mse
#' @export lm.mse
#' @param X covariates (n times p matrix, n: number of entries, p: number of covariates)
#' @param Y response (vector with n entries)
#' @param mod fitted model from lm.cv or csuv. Only provide mod or est.b
#' @param est.b estimated coefficient (with intercept). Only provide mod or est.b
#' @param log.level log level to set. Default is NULL, which means no change in log level. See the function CSUV::set.log.level for more details
#' @return the value of estimated mean square error
#' @examples
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' compare.mod = lm.compare.method(X, Y, intercept = FALSE)
#' lm.mse(X, Y, est.b = compare.mod)
lm.mse <- function(X, Y, mod = NULL, est.b = NULL, log.level = NULL) {
  if (!is.null(log.level)){
    set.log.level(log.level)
  }
  pred <- lm.predict(X, mod = mod, est.b = est.b)
  if (is.null(mod) == is.null(est.b)) {
    stop("wrong param for function lm.mse")
  }
  return(colMeans((Y - pred)^2))
}

lm.predict <- function(X, mod = NULL, est.b = NULL) {
  # est.b is with length p+1 (first one is intercept)

  # check param
  if (is.null(mod) == is.null(est.b)) {
    stop("wrong param for function lm.predict")
  }

  # set beta
  if (is.null(est.b)) {
    est.b <- mod$est.b
  }
  if (is.vector(est.b)) {
    return(X %*% est.b[-1] + est.b[1])  # return prediction for only one model
  } else{
    est.b <- as.matrix(est.b)
    b0 <- matrix(est.b[, 1], nrow = nrow(X), ncol = nrow(est.b), byrow = T)
    return(X %*% t(est.b[, -1, drop = FALSE]) + b0) # for several models
  }

}

## ==== helper ====
# there is some problem with cvrelaxo so we use cvrelaxo1 here
cvrelaxo1 <- function(X, Y, K = 5, phi = seq(0, 1, length = 10), max.steps = min(2 * length(Y), 2 * ncol(X)), fast = TRUE, keep.data = TRUE, warn = TRUE) {
  Y <- as.numeric(Y)
  if (warn) {
    if (abs(mean(Y)) > 0.01 * stats::sd(Y))
      warning("response variable not centered")
    if (any(abs(apply(X, 2, mean)) > 0.01 * apply(X, 2, stats::sd)))
      warning("predictor variables not centered")
    if (stats::sd(as.numeric(apply(X, 2, stats::sd))) > 0.001)
      warning("predictor variables not scaled")
  }
  n <- length(Y)
  index <- sample(rep(1:K, each = ceiling(n / K)), n, replace = FALSE)
  losscv <- rep(0, length = length(phi) * (max.steps - 1))
  for (k in 1:K) {
    rel <- relaxo::relaxo(X[index != k, ], Y[index != k], phi = phi,
                  fast = fast, keep.data = FALSE, warn = FALSE,
                  max.steps = max.steps)
    pred <- X[index == k, ] %*% t(rel$beta)
    losscv[seq_len(ncol(pred))] <- losscv[seq_len(ncol(pred))] + apply(sweep(pred, 1, Y[index == k])^2, 2, mean) / K
    if (length(losscv) > ncol(pred))
      losscv[(ncol(pred) + 1):length(losscv)] <- Inf
  }
  relall <- relaxo::relaxo(X, Y, phi = phi, fast = fast,
                   keep.data = keep.data, warn = FALSE)
  select <- which.min(losscv[seq_len(nrow(relall$beta))]) # quick fix...
  relall$beta <- relall$beta[select, , drop = FALSE]
  relall$lambda <- relall$lambda[select]
  relall$phi <- relall$phi[select]
  return(relall)
}
