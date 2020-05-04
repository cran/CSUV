cv.for.one.method <- function(X, Y, intercept, flds, fit.method, is.keep.all, progress.percent = NULL) {
  # helper function of "cv.for.methods", which creates folds of cv fit for one method
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: whether the model should have intercept
  #   fit.method: fitting method
  #   flds: test folds
  #   is.keep.all: whether we keep the fitted models, when ols fit is not available
  #
  # Returns:
  #   a list (index by folds) of list (estimated beta and mse)

  futile.logger::flog.info("start fitting cv with method %s", attr(fit.method, "name"))

  result <- list()
  if (!is.null(csuv.env$is.shiny) && csuv.env$is.shiny && !is.null(progress.percent)) {
    result <- lapply(flds, function(fld) {
      cv.result <- cv.helper.one.fld.one.method.with.check(X, Y, intercept, fld, fit.method, is.keep.all)
      shiny::incProgress(progress.percent / length(flds))
      return(cv.result)
    })
  } else{
    result <- lapply(flds, function(fld)
      cv.helper.one.fld.one.method.with.check(X, Y, intercept, fld, fit.method, is.keep.all))
  }
  futile.logger::flog.info("finished fitting cv with method %s", attr(fit.method, "name"))

  return(result)
}

## ------ cv helpers ------
cv.helper.one.fld.one.model <- function(X, Y, intercept, fld, est.b = NULL,
                                      fit.method = NULL, with.ols = TRUE) {
  # helper function which creates one fold cv fit for one specific model or method
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: whether the model should have intercept
  #   fld: test fold
  #   est.b: estimated beta. If it is not NULL, then the fit is for a specific fitted model
  #   fit.method: fitting method. If it is not NULL, then the fit is for a specific method
  #
  # Returns:
  #   a list (ols.cv, cv) of list (estimated beta and mse)

  is.cv.on.fitted.model <- (!is.null(est.b) && is.null(fit.method) && with.ols)
  is.cv.on.method <- (is.null(est.b) && !is.null(fit.method))

  cv.fit.result <- list()
  if (is.cv.on.method) {
    cv.fit.result$cv <- fit.method(X[-fld, , drop = FALSE], Y[-fld], intercept)
    est.b <- cv.fit.result$cv$est.b

    if (with.ols) {
      if (is.vector(est.b)) {
        if (length(sum(est.b != 0)) < length(Y[-fld])) {
          cv.fit.result$ols.cv <- lm.ols.refit(X[-fld, , drop = FALSE], Y[-fld],
                                              intercept, est.b)
        }
      } else { # if est.beta is a matrix
        sel.i <- which(rowSums(est.b != 0) < length(Y[-fld])) # make sure there are enough sample to refit
        cv.fit.result$ols.cv <- lm.ols.refit(X[-fld, , drop = FALSE], Y[-fld],
                                             intercept, est.b[sel.i, , drop = FALSE])
      }
    }
  } else if (is.cv.on.fitted.model) {
    cv.fit.result$ols.cv <- lm.ols.refit(X[-fld, , drop = FALSE], Y[-fld], intercept, est.b)
  } else{
    futile.logger::flog.error("cv.for.one.fld.one.model: wrong argument(s)")
    stop("wrong argument for cv.for.one.fld.one.model")
  }

  # get result with mse
  result <- list()
  for (cv.name in names(cv.fit.result)) {
    result[[cv.name]] <- list(est.b = cv.fit.result[[cv.name]]$est.b,
                              mse = lm.mse(X[fld, , drop = FALSE], Y[fld],
                                           mod = cv.fit.result[[cv.name]]))
  }
  return(result)
}

cv.helper.one.fld.one.method.with.check <- function(X, Y, intercept, fld,
                                                  fit.method, is.keep.all) {
  # helper function of "cv.for.one.cv.method", which creates one fold cv fit for one specific method with check
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: whether the model should have intercept
  #   fld: test fold
  #   fit.method: fitting method
  #   is.keep.all: whether we keep the fitted models, when ols fit is not available
  #
  # Returns:
  #   a list (estimated beta and mse)
  cv.result <- cv.helper.one.fld.one.model(X, Y, intercept = intercept, fld = fld,
                                           fit.method = fit.method, with.ols = TRUE)

  b <- cv.result$cv$est.b
  e <- cv.result$cv$mse
  ols.b <- NULL
  ols.e <- NULL

  if (!is.null(cv.result$ols.cv)) { # if can fit ols, use it
    ols.b <- cv.result$ols.cv$est.b
    ols.e <- cv.result$ols.cv$mse
  } else if (is.keep.all) { # if cannot fit ols and we need to keep all methods, use the raw one for ols ones
    ols.b <- b
    ols.e <- e
  }
  if (any(is.na(b)) || any(is.na(ols.b))) {
    stop("estimated beta cannot be NA")
  }
  return(list(est.b = b, mse = e, ols.est.b = ols.b, ols.mse = ols.e))
}

## ======= create folds =========
get.flds <- function(X, fit.percent = NULL, num.repeat = NULL,
                     is.k.fld.cv = FALSE, num.fld = NULL,
                     min.train.size = 0, current.flds = NULL) {
  # create a list of test folds, with guarentees that each covariate in each train fold is not constant,
  # and the train fold size is not smaller than min.train.size
  #
  # Args:
  #   X: covariates
  #   fit.percent: percentage of obs used in fitting
  #   num.repeat: number of time to get samples
  #   min.train.size: minimum size of each train fold
  #   current.flds: current folds
  #
  # Returns:
  #   a list of test folds

  if (is.k.fld.cv) {
    return(get.good.flds(X, num.fld))
  } else {
    n <- nrow(X)
    test.size <- min(n - min.train.size, round((1 - fit.percent) * n))

    current.num.repeat <- length(current.flds)
    num.from.current <- if (length(current.flds) == 0) { NULL } else { 1:min(current.num.repeat, num.repeat)}
    num.to.create <- if (num.repeat <= current.num.repeat) { NULL } else { (current.num.repeat + 1):num.repeat}

    return(c(lapply(num.from.current, function(i) current.flds[[i]]),
              lapply(num.to.create, function(i) get.good.samples(X, test.size, seed = NULL)))) #i
  }
}

get.flds.for.methods <- function(flds, current.fld.length.for.method) {
  # create a list of test folds for a method
  #
  # Args:
  #   flds: folds for all methods
  #   current.fld.length.for.method: length of folds that already fitted
  #
  # Returns:
  #   a list of test folds with length = max(0, current.fld.length.for.method, length(flds))

  target.fld.length <- length(flds)
  if (target.fld.length <= current.fld.length.for.method) {
    return()
  } else{
    return(lapply((current.fld.length.for.method + 1):target.fld.length,
                   function(i) flds[[i]]))
  }
}

## ------ create fold helpers ------
get.good.samples <- function(X, test.size, seed = NULL) {
  # helper function of "create.flds", which creates delete-k cv samples with guarentee that each covariate in each train fold is varied
  #
  # Args:
  #   X: covariates
  #   test.size: size of each test fold (k)
  #
  # Returns:
  #   a vector of test folds

  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- ncol(X)
  varied.i <- which(sapply(1:p, function(i) {
    if (is.numeric(X[, i])) {
      return(max(X[, i]) != min(X[, i]))
    } else if (is.factor(X[, i])) {
      return(nlevels(X[, i]) > 1)
    }

  }))
  counter <- 20
  repeat{
    n <- nrow(X)
    fld <- sample(1:n, test.size)
    if (!counter || !any(sapply(varied.i, function(i) {
      if (is.numeric(X[-fld, i])) {
        return(max(X[-fld, i]) == min(X[-fld, i]))
      } else{
        return(nlevels(X[-fld,i]) > 1)
      }

      }))) {
      return(fld)
    }
    counter <- counter - 1
    futile.logger::flog.info("get.good.samples: redraw")
  }
}

get.good.flds <- function(X, num.fld) {
  # helper function of "create.flds", which creates k-fold CV samples with guarentee
  # that each covariate in each train fold is varied
  #
  # Args:
  #   X: covariates
  #   num.fld: number of flds
  #
  # Returns:
  #   a vector of test folds

  varied.i <- which(sapply(seq_len(ncol(X)), function(i) max(X[, i]) != min(X[, i])))
  n <- nrow(X)
  repeat{
    flds <- caret::createFolds(1:n, k = num.fld, list = TRUE)
    if (!any(sapply(flds, function(fld)
      any(sapply(varied.i, function(i) max(X[-fld, i]) == min(X[-fld, i])))))) {
      return(flds)
    }
    futile.logger::flog.info("get.good.samples: redraw")
  }
}

cv.for.methods <- function(X, Y, intercept,
                         fit.percent, num.repeat,
                         fit.methods, is.keep.all,
                         current.fit = NULL, num.core = 1) {
  # the function creates delete k cv fit for methods
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: whether the model should have intercept
  #   fit.percent: percentage of obs used in fitting
  #   num.repeat: number of iterations
  #   fit.methods: fitting methods
  #   is.keep.all: whether we keep the fitted models, when ols fit is not available
  #   current.fit: current fits
  #
  # Returns:
  #   a list of: 1. list (index by methods) of list (index by folds) of list (estimated beta and mse)
  #              2. flds

  cv.name <- paste0(round(fit.percent, 2), ":", round(1 - fit.percent, 2),
                    " cv for ", num.repeat, " times")
  futile.logger::flog.info("start fitting %s", cv.name)

  current.flds <- NULL
  if (!is.null(current.fit)) {
    current.flds <- current.fit$flds
  }
  flds <- get.flds(X, fit.percent = fit.percent, num.repeat = num.repeat,
                  min.train.size = 0, current.flds = current.flds)

  result <- NULL
  if (num.core > 1) {
    cl <- parallel::makeCluster(min(parallel::detectCores() - 1, length(fit.methods), num.core),
                      outfile = "parallel_log.txt")
    parallel::clusterExport(cl = cl, list("get.flds.for.methods", "cv.for.one.method"),
                  envir = environment())
    result <- parallel::parLapply(cl, names(fit.methods), function(fit.method.name, param) {
      method.flds <- get.flds.for.methods(param$flds, current.fld.length.for.method = length(param$current.fit$unique.mod[[fit.method.name]]))
      return(cv.for.one.method(param$X, param$Y, param$intercept, method.flds, param$fit.methods[[fit.method.name]], param$is.keep.all, progress.percent = NULL))
    }, param = list(X = X, Y = Y, intercept = intercept,
                    flds = flds, current.fit = current.fit,
                    fit.methods = fit.methods, is.keep.all = is.keep.all))
    parallel::stopCluster(cl)
  } else{
    result <- lapply(names(fit.methods), function(fit.method.name) {
      if (!is.null(csuv.env$is.shiny) && csuv.env$is.shiny) {
        shiny::incProgress(amount = 0, detail = paste("\n fitting", fit.method.name))
      }
      method.flds <- get.flds.for.methods(flds, current.fld.length.for.method = length(current.fit$unique.mod[[fit.method.name]]))
      return(cv.for.one.method(X, Y, intercept, method.flds,
                               fit.methods[[fit.method.name]], is.keep.all,
                               progress.percent = 1 / length(fit.methods)))
    })
  }

  names(result) <- names(fit.methods)

  futile.logger::flog.info("finished fitting %s", cv.name)
  return(list(cv.result = result,
              flds = flds))
}
