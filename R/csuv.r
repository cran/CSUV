#################################################
## This file includes the main functions for CSUV
#################################################
csuv.env <- new.env()

#' Print the coefficients of csuv
#' @export
#' @param x output of csuv()
#' @param ... additional arguments to "print"
#' @return return value from print(x$est.b)
print.csuv <- function(x, ...) {
  print(x$est.b, ...)
}

#' Get the fitted results from Combined Selection and Uncertainty Visualiser (CSUV) method
#' @export csuv
#' @import parallel
#' @import doParallel
#' @import relaxo
#' @import datasets
#' @importFrom grDevices rainbow rgb dev.list dev.off
#' @importFrom graphics abline axis boxplot legend par plot points rect text plot.new
#' @importFrom stats complete.cases median sd lm coef quantile
#' @importFrom utils capture.output head read.csv
#' @param X covariates (n times p matrix, n: number of entries, p: number of covariates)
#' @param Y response (vector with n entries)
#' @param intercept TRUE to fit the data with an intercept, FALSE to fit the data without an intercept
#' @param method.names vector of method names to be used in CSUV. Choose among "lasso", "elastic", "relaxo", "mcp" and "scad". Default is to use all methods listed above
#' @param coef.est.method method to estimate the coefficients of covariates after variable selection. User can provide his/her function. Default is ordinary least square
#' @param B number of subsampling. Default is 100
#' @param fit.percent percentage of observations used in fitting in CSUV
#' @param q percentile of fitted models used per each subsampling in CSUV, according to the selection criterion on out-of-sample data in ascending order. Default is q = 0 (only the fitted model with the lowest MSE in a subsampling data is used)
#' @param selection.criterion = c("mse", "ebic"). Measure to select fitted models in subsampling dataset. "mse" is mean square error and "ebic" is extended BIC. Default is mse
#' @param num.core number of cores to use. Default is 1 (i.e. no parallel running)
#' @param all.fits (optional) all fitted models. If all.fits is provided, then CSUV will use the fitted models in all.fitted instead of fitting using subsampling data
#' @param log.level log level to set. Default is NULL, which means no change in log level. See the function CSUV::set.log.level for more details
#' @return a list, which includes estimated coefficients (est.b), subsampling fitted models (mod.collection), number of times a method is selected (method.freq), relative frequency of each covariate (variable.freq), covariates ordered by relative frequency (variable.order).
#' @examples
#' \donttest{
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' mod.0 = csuv(X, Y, intercept = FALSE, q = 0, method.names = NULL)
#' print(mod.0)
#' mod.5 = csuv(X, Y, intercept = FALSE, q = 5, all.fits = mod.0$all.fits)
#' print(mod.5)
#' }
csuv <- function(X, Y, intercept, method.names = NULL,
                 coef.est.method = lm.ols,
                 B = 100, q = 0,
                 fit.percent = 0.5,
                 selection.criterion = "mse",
                 num.core = 1,
                 all.fits = NULL,
                 log.level = NULL) {
  # wrapper function of the main algorithm
  if (!is.null(log.level)){
    set.log.level(log.level)
  }
  if (is.null(method.names)) {
    method.names <- names(get.fit.methods())
  }
  mod <- get.csuv.mod(X = X, Y = Y, intercept = intercept,
                      method.names = method.names,
                      coef.est.method = coef.est.method,
                      B = B, q = q, fit.percent = fit.percent,
                      selection.criterion = selection.criterion,
                      num.core = num.core,
                      all.fits = all.fits)
  class(mod) <- "csuv"
  return(mod)
}


get.csuv.mod <- function(X, Y, intercept,
                         method.names, coef.est.method,
                         B, q, fit.percent,
                         selection.criterion,
                         num.core,
                         all.fits = NULL) {
  # the main function of csuv
  #
  #   X: covariates
  #   Y: response
  #   intercept: intercept
  #   method.names: names of methods to use
  #   coef.est.method: refitting method
  #   B: number of iteration in CV
  #   q: percentile of fitted models from cv to use
  #
  # Returns:
  #   a list (variable selection frequency, solution path)
  if (is.null(all.fits)) {
    all.fits <- get.csuv.unique.fit(X = X, Y = Y, intercept = intercept,
                                    method.names = method.names,
                                    B = B,
                                    fit.percent = fit.percent,
                                    current.fit = NULL,
                                    num.core = num.core)
  }
  csuv.result <- get.csuv.final.mod(X = X, Y = Y, intercept = intercept,
                                    unique.fit = all.fits,
                                    selection.criterion = selection.criterion,
                                    coef.est.method = coef.est.method,
                                    q = q,
                                    method.names = method.names,
                                    B = B)
  csuv.result$all.fits <- all.fits
  csuv.result$param <- list(selection.criterion = selection.criterion,
                            coef.est.method = coef.est.method,
                            fit.percent = fit.percent,
                            q = q,
                            method.names = method.names,
                            B = B)
  return(csuv.result)
}

#' Helper function, please do not use it
#' @export get.csuv.unique.fit
#' @param X covariates (n times p matrix, n: number of entries, p: number of covariates)
#' @param Y response (vector with n entries)
#' @param intercept TRUE to fit the data with an intercept, FALSE to fit the data without an intercept
#' @param method.names vector of method names to be used in CSUV. Choose among "lasso", "elastic", "relaxo", "mcp" and "scad". Default is to use all methods listed above
#' @param B number of subsampling. Default is 100
#' @param fit.percent percentage of observations used in fitting in CSUV
#' @param num.core number of cores to use. Default is 1 (i.e. no parallel running)
#' @param current.fit (optional) all fitted models
#' @return a list of current fit
get.csuv.unique.fit <- function(X, Y, intercept, method.names, B,
                                fit.percent,
                                current.fit = NULL,
                                num.core = 1) {
  # the main function of the new method
  #
  #   X: covariates
  #   Y: response
  #   intercept: intercept
  #   method.names: names of methods to use
  #   B: number of iteration in CV
  #
  # Returns:
  #   a list (methods) of list (flds) of list (unique models (with coefficient), mse)
  path.fit.methods <- get.path.fit.methods()
  cv.fit <- cv.for.methods(X, Y, intercept,
                           fit.percent, num.repeat = B,
                           fit.methods = path.fit.methods[method.names],
                           is.keep.all = FALSE, current.fit = current.fit,
                           num.core = num.core)
  unique.mod <- unique.mod.for.methods(cv.fit$cv.result, is.ols = TRUE)

  current.method.names <- union(names(current.fit$unique.mod), names(unique.mod))
  current.mod <- lapply(current.method.names,
                        function(method.name) c(current.fit$unique.mod[[method.name]], unique.mod[[method.name]]))

  names(current.mod) <- current.method.names

  current.fit <- list(unique.mod = current.mod,
                      flds = if (length(cv.fit$flds) > length(current.fit$flds)) cv.fit$flds else current.fit$flds)

  attr(current.fit, "type") <- "csuv.all.fits"
  return(current.fit)
}

is.csuv.all.fits <- function(obj) {
  return(!is.null(attr(obj, "type")) && attr(obj, "type") == "csuv.all.fits")
}

is.csuv.fit <- function(obj) {
  return(class(obj) == "csuv")
}

#' Helper function, please do not use it
#' @export get.csuv.final.mod
#' @param X covariates (n times p matrix, n: number of entries, p: number of covariates)
#' @param Y response (vector with n entries)
#' @param intercept TRUE to fit the data with an intercept, FALSE to fit the data without an intercept
#' @param unique.fit from get.csuv.unique.fit
#' @param selection.criterion = c("mse", "ebic"). Measure to select fitted models in subsampling dataset. "mse" is mean square error and "ebic" is extended BIC. Default is mse
#' @param coef.est.method method to estimate the coefficients of covariates after variable selection. User can provide his/her function. Default is ordinary least square
#' @param q percentile of fitted models used per each subsampling in CSUV, according to the selection criterion on out-of-sample data in ascending order. Default is q = 0 (only the fitted model with the lowest MSE in a subsampling data is used)
#' @param method.names vector of method names to be used in CSUV. Choose among "lasso", "elastic", "relaxo", "mcp" and "scad". Default is to use all methods listed above
#' @param B number of subsampling. Default is 100
#' @return a list of current fit
get.csuv.final.mod <- function(X, Y, intercept, unique.fit,
                               selection.criterion,
                               coef.est.method = lm.ols,
                               q, method.names, B) {
  # the main function of the new method
  #
  #   X: covariates
  #   Y: response
  #   intercept: intercept
  #   unique.fit: unique fit from "get.csuv.unique.fit"
  #   coef.est.method: refitting method
  #   q: percentile of fitted models from cv to use
  #
  # Returns:
  # a list (methods) of list (flds) of list (unique models (with coefficient), mse)
  # filter the models and calculate the variable selection frequency
  sel.mods.n.sel.freq <- get.best.mod(unique.fit,
                                      selection.criterion = selection.criterion,
                                      is.ols = TRUE, q = q,
                                      method.names = method.names,
                                      B = B)
  var.sel.freq.n.order <- csuv.var.sel.freq.n.sol.path(sel.mods.n.sel.freq)

  # get the final models
  csuv.s <- csuv.thr.by.size(X, Y, intercept,
                             var.order = var.sel.freq.n.order$var.order,
                             size.thr = stats::median(rowSums(sel.mods.n.sel.freq$sel.mod[, -1, drop = FALSE] != 0)),
                             coef.est.method)

  csuv.m <- csuv.thr.by.freq(X, Y, intercept,
                             var.sel.freq = var.sel.freq.n.order$var.sel.freq,
                             freq.thr = 0.5,
                             coef.est.method)
  csuv.ebic <- csuv.by.ebic(X, Y, intercept,
                            var.order = var.sel.freq.n.order$var.order,
                            size.thr = sum(var.sel.freq.n.order$var.sel.freq >= 0.1))
  est.b <- rbind(csuv.m = csuv.m,
                csuv.s = csuv.s,
                csuv.ebic = csuv.ebic)

  colnames(est.b) <- colnames(sel.mods.n.sel.freq$sel.mod)
  colnames(est.b)[1] <- "intercept"
  return(list(est.b = est.b,
               mod.collection = sel.mods.n.sel.freq$sel.mod,
               method.freq = sel.mods.n.sel.freq$sel.method.freq,
               variable.freq = var.sel.freq.n.order$var.sel.freq,
               variable.order = var.sel.freq.n.order$var.order))
}
