#################################################
## This file includes the helper functions for CSUV
#################################################

## ====== filtering fitted models =======
get.best.mod <- function(cv.results, selection.criterion, is.ols,
                       q, method.names, B) {
  # a function to get the best k (size) fitted models in folds
  #
  # Args:
  #   cv.results: cv results from "cv.for.methods"
  #   is.ols: whether ols or original cv results should be used
  #   fld.indexes: fld to use. If is.null, then use all
  #
  # Returns:
  #   a list (index by folds) of list (mods, selection freq)

  # change from list (method) of list (fold) of list (est.b/mse) of
  # matrix (path*p+1) to list (fold) of matrix (path*p+1)
  est.b <- reshape.cv.result(cv.results$unique.mod,
                             item.name = ifelse(is.ols, "ols.est.b", "est.b"),
                             method.names, B)
  mse <- reshape.cv.result(cv.results$unique.mod,
                           item.name = ifelse(is.ols, "ols.mse", "mse"),
                           method.names, B)

  # get ebic / bic / mse
  criterion.val <- mse
  if (selection.criterion == "ebic") {
    criterion.val <- lapply(seq_len(length(mse)), function(fld.idx) {
      ebic.val <- get.ebic.value.from.mse(mse = mse[[fld.idx]],
                                         n = length(cv.results$flds[[fld.idx]]),
                                         est.b = est.b[[fld.idx]])
      return(ebic.val)
    })

  } else if (selection.criterion == "bic") {
    criterion.val <- lapply(seq_len(length(mse)), function(fld.idx) {
      bic.val <- get.bic.value.from.mse(mse = mse[[fld.idx]],
                                        n = length(cv.results$flds[[fld.idx]]),
                                        est.b = est.b[[fld.idx]])
      return(bic.val)
    }
    )
  }

  # calculate the number of fitted models to retain
  num.retain <- 1
  if (q != 0) {
    avg.mod.size <- mean(sapply(mse, length))
    num.retain <- round(q / 100 * avg.mod.size)
  }

  ## ===== get the fitted models with lowest mse/bic ========
  result <- lapply(1:B, function(fld.idx)
    get.best.mod.helper.one.fld(method.names = method.names,
                                est.b[[fld.idx]],
                                criterion.val[[fld.idx]],
                                num.retain))
  return(list(sel.mod = do.call(rbind,
                                 lapply(result, function(x) x$sel.mod)),
               sel.method.freq <- do.call(rbind,
                                         lapply(result, function(x)
                                           x$sel.method.freq))))
}

## ---- helper functions to get the best fitted models -----
reshape.cv.result <- function(unique.mod, item.name, method.names, B) {
  # helper functions which re-structure the unique.mod
  #
  # Args:
  #   unique.mod: cv results from "cv.for.methods"
  #   item.name: item to get (e.g. est.b)
  #   fld.indexes: fld to use
  #
  # Returns:
  #   a list (index by folds) of matrix (index by methods) of the specific item

  if (is.vector(unique.mod[[1]][[1]][[item.name]])) {
    result <- lapply(1:B,
                    function(fld.idx) do.call(c, lapply(method.names,
                                                        function(method.name)
                                                          unique.mod[[method.name]][[fld.idx]][[item.name]])))
    return(result)
  } else{
    return(lapply(1:B,
                   function(fld.idx)
                     do.call(rbind, lapply(method.names,
                                           function(method.name) {
                                             x <- unique.mod[[method.name]][[fld.idx]][[item.name]]
                                             rownames(x) <- rep(method.name, nrow(x))
                                             return(x)
                                             }))))
  }

}

get.best.mod.helper.one.fld <- function(method.names, est.b, mse, num.retain) {
  # a helper function for "get.best.mod.given.size" to get the
  # best k (num.retain) fitted models in one fold
  #
  # Args:
  #   method.names: method names
  #   est.b: estimated coefficients of the fitted models
  #   mse: mse
  #   num.retain: number of fitted models to retain
  #
  # Returns:
  #   a list (mods, selection freq)

  num.retain <- min(length(mse), num.retain) # make sure num of retain is not larger than number of models
  e.order <- order(mse, rowSums(est.b != 0))
  sel.mod <- est.b[e.order[1:num.retain], , drop = FALSE]

  mse.thr <- mse[e.order[num.retain]]
  candidate.b <- est.b[which(mse <= mse.thr), , drop = FALSE]
  candidate.e <- mse[which(mse <= mse.thr)]
  last.size <- sum(est.b[e.order[num.retain], ] != 0)
  sel.method.names <- rownames(candidate.b[which(candidate.e < mse.thr | rowSums(candidate.b != 0) == last.size), , drop = FALSE])
  sel.method.freq <- sapply(method.names, function(method.name)
    sum(sel.method.names == method.name))

  return(list(sel.mod = sel.mod, sel.method.freq = sel.method.freq))
}

## ======== calculate variable selection frequency and get solution path =============
csuv.var.sel.freq.n.sol.path <- function(sel.mods.n.sel.freq) {
  # a function to calculate variable selection frequency and get solution path
  #
  # Args:
  #   sel.mod from "get.best.mod"
  #
  # Returns:
  #   a list (variable selection frequency, soluation path)
  sel.mods <- sel.mods.n.sel.freq$sel.mod
  # sel.freq <- sel.mods.n.sel.freq$sel.method.freq
  var.sel.freq <- sapply(2:ncol(sel.mods), # exclude intercept
                         function(x) max(mean(sel.mods[, x] > 0),
                                          mean(sel.mods[, x] < 0)))
  mean.beta <- abs(colMeans(sel.mods))[-1]# exclude intercept
  var.order <- order(-var.sel.freq, -mean.beta)
  return(list(var.sel.freq = var.sel.freq,
               var.order = var.order))
}

## ====== get the final fitted models after filtering ==========

csuv.thr.by.size <- function(X, Y, intercept, var.order, size.thr, coef.est.method) {
  # a function to get the final fitted model using size threshold
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: intercept
  #   var.order: variable ordering
  #   size.thr: number of variable to retain
  #   coef.est.method: refitting method
  #
  # Returns:
  #   a fitted model

  est.i <- if (size.thr > 0) { est.i <- var.order[1:size.thr] } else c()
  return(csuv.refit(X, Y, intercept, est.i, coef.est.method))
}

csuv.thr.by.freq <- function(X, Y, intercept, var.sel.freq, freq.thr, coef.est.method) {
  # a function to get the final fitted model using frequency threshold
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: intercept
  #   var.order: variable ordering
  #   freq.thr: minimum frequency in order to retain a variable
  #   coef.est.method: refitting method
  #
  # Returns:
  #   a fitted model

  est.i <- which(var.sel.freq >= freq.thr)
  return(csuv.refit(X, Y, intercept, est.i, coef.est.method))
}

csuv.by.ebic <- function(X, Y, intercept, var.order, size.thr) {
  # a function to get the final fitted model using size threshold
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: intercept
  #   var.order: variable ordering
  #   size.thr: number of variable to retain
  #   coef.est.method: refitting method
  #
  # Returns:
  #   a fitted model

  p <- ncol(X)
  est.b <- do.call(rbind, lapply(1:size.thr, function(size) {
    b <- rep(0, p)
    b[var.order[1:size]] <- 1
    return(b)
  }))
  est.b <- cbind(intercept, est.b)
  est.b <- lm.ols.refit(X, Y, intercept, est.b)$est.b
  mse <- lm.mse(X, Y, est.b = est.b)

  n <- length(Y)
  ebic.val <- get.ebic.value.from.mse(mse, n, est.b)
  return(est.b[which.min(ebic.val), ])
}

csuv.refit <- function(X, Y, intercept, est.i, coef.est.method) {
  # a helper function for "csuv.thr.by.freq" and "csuv.thr.by.size"
  # to refit the selected variables
  #
  # Args:
  #   X: covariates
  #   Y: response
  #   intercept: intercept
  #   est.i: selected covariate indexes
  #   coef.est.method: refitting method
  #
  # Returns:
  #   a fitted model

  p <- ncol(X)
  new.b <- rep(0, (p + 1))
  if (length(est.i)) {
    new.b[c(0, est.i) + 1] <- coef.est.method(X[, est.i, drop = FALSE], Y,
                                              intercept = intercept)$est.b
  }
  return(new.b)
}

## ====== sort and remove duplication ==========
unique.mod.for.methods <- function(cv.methods.result, is.ols) {
  # a function which removes duplicated models (in terms of selection)
  #
  # Args:
  #   cv.folds.result: cv result of methods with folds
  #   is.ols: whether the ols or the original fitting result should be used
  #   fld.idx: which fld index to use. If omitted, all results will be used
  #
  # Returns:
  #   a list (methods) of list (flds) of list(est.b, mse)

  return(lapply(cv.methods.result, function(cv.folds.result)
    unique.mod.helper.one.method(cv.folds.result, is.ols)))
}

unique.mod.helper.one.method <- function(cv.flds.result, is.ols) {
  # a helper function of "unique.mod.for.methods" which removes duplicated models (in terms of selection)
  #
  # Args:
  #   cv.folds.result: cv result of a method with folds
  #   is.ols: whether the ols or the original fitting result should be used
  #   fld.idx: which fld index to use. If omitted, all results will be used
  #
  # Returns:
  #   a list (flds) of list(est.b, mse)

  est.b.name <- ifelse(is.ols, "ols.est.b", "est.b")
  mse.name <- ifelse(is.ols, "ols.mse", "mse")

  return(lapply(cv.flds.result, function(cv.fld.result) {
    est.b <- cv.fld.result[[est.b.name]]
    attr(est.b, "name") <- est.b.name
    mse <- cv.fld.result[[mse.name]]
    attr(mse, "name") <- mse.name
    get.unique.mod(est.b = est.b, mse = mse, is.ols = is.ols)
    }))
}

get.unique.mod <- function(est.b, mse, is.ols) {
  # helper function of "unique.mod.helper.one.method", which removes duplicated models (in terms of selection)
  #
  # Args:
  #   est.b: estimated coefficients
  #   mse: mse
  #
  # Returns:
  #   a list (est.b, mse)

  # if not ols, sort the est.b so that the one with lowest mse will be selected
  if (!is.ols) {
    e.order <- order(mse, rowSums(est.b != 0)) # order by mse, then by model size
    est.b <- est.b[e.order, ]
    mse <- mse[e.order, ]
  }

  if(is.null(dim(est.b))){
    return (list(est.b = matrix(est.b, nrow = 1),
                 mse = mse))
  }
  non.duplicate.i <- which(!duplicated(est.b != 0))
  result <- list(est.b = est.b[non.duplicate.i, , drop = FALSE],
                 mse = mse[non.duplicate.i])
  if (any(is.na(result$est.b)) || any(is.na(result$mse))) {
    stop("get.unique.mod: estimated beta or mse cannot be NA")
  }
  names(result) <- c(attr(est.b, "name"), attr(mse, "name"))

  return(result)
}

get.ebic.value.from.mse <- function(mse, n, est.b) {
  p <- ncol(est.b[, -1, drop = FALSE]) # number of covariate, first one is intercept
  k <- rowSums(est.b[, -1, drop = FALSE] != 0) # selected model size
  return(get.bic.value.from.mse(mse, n, est.b) + log(choose(p, k)))
}

get.bic.value.from.mse <- function(mse, n, est.b) {
  k <- rowSums(est.b[, -1, drop = FALSE] != 0) # selected model size
  bic.value <- n * log(mse) + log(n) * k
  return(bic.value)
}
