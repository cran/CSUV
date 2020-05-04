#' Confidence interval-like interval for uncertainty illustration
#' @export csuv.ci
#' @param csuv.fit fitted results from CSUV::csuv()
#' @param level significance level
#' @param type the type of the interval. When the type is "original", all estimated coefficients are used to calculate the interval. When the type is "conditional", only non-zero estimated coefficients are used. The type "conditional.1" is still in experimental stage, please do not use. Default is "original"
#' @param log.level log level to set. Default is NULL, which means no change in log level. See the function CSUV::set.log.level for more details
#' @return a matrix. Each column represents an interval for a corresponding covariate
#' @examples
#' \donttest{
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' mod.0 = csuv(X, Y, intercept = FALSE, q = 0, method.names = NULL)
#' print(csuv.ci(mod.0, level = 0.1, type = "original"))
#' }
csuv.ci <- function(csuv.fit, level, type = "original", log.level = NULL) {
  if (!is.null(log.level)){
    set.log.level(log.level)
  }
  mods <- NULL
  if (type == "original") {
    mods <- csuv.fit$mod.collection[, -1] # first is intercept
  } else if (type == "conditional") {
    mods <- csuv.fit$mod.collection[, -1]
    mods[which(mods == 0)] <- NA
  } else if (type == "conditional.1") {
    mods <- csuv.fit$mod.collection[, -1]
    pos.i <- which(colMeans(mods > 0) > colMeans(mods < 0))
    neg.i <- which(colMeans(mods > 0) < colMeans(mods < 0))
    if (length(pos.i)) mods[, pos.i][which(mods[, pos.i] < 0)] <- NA
    if (length(neg.i)) mods[, neg.i][which(mods[, neg.i] > 0)] <- NA
    mods[which(mods == 0)] <- NA
  }
  p <- ncol(mods)

  li <- NULL
  ui <- NULL
  if (type == "conditional") {

    li <- sapply(1:p, function(i) {
      stats::quantile(mods[which(mods[, i] != 0), i], level / 2,
                      na.rm = TRUE)
    })
    ui <- sapply(1:p, function(i) {
      stats::quantile(mods[which(mods[, i] != 0), i], 1 - level / 2,
                      na.rm = TRUE)
    })
  } else{
    li <- sapply(1:p, function(i)
      stats::quantile(mods[, i], level / 2, na.rm = TRUE))
    ui <- sapply(1:p, function(i)
      stats::quantile(mods[, i], 1 - level / 2, na.rm = TRUE))
  }

  names(li) <- colnames(mods)
  li[which(is.na(li))] <- 0
  ui[which(is.na(ui))] <- 0
  return(rbind(li = li,
                ui = ui))
}

### ==== LPR helper ====
lm.PR <- function(x, y, intercept, est.b, ...) {
  PR.mod <- HDCI::PartRidge(x = x, y = y, intercept = intercept,
                            varset = est.b[-1] != 0,
                            lambda2 = 1 / length(y), ...)
  return(c(PR.mod$beta0, PR.mod$beta))
}

get.PR.mod.collection <- function(x, y, intercept, flds,
                                  csuv.mod.collection,
                                  variable.freq, num.core = 1) {
  result <- NULL
  if (num.core > 1) {
    cl <- parallel::makeCluster(min(parallel::detectCores() - 1, num.core),
                      outfile = "parallel_log.txt")
    parallel::clusterExport(cl = cl, list("lm.PR"), envir = environment())
    result <- parallel::parLapply(cl, seq_len(length(flds)), function(fld.i, param) {
      fld <- param$flds[[fld.i]]
      return(lm.PR(x = param$x[-fld, ], y = param$y[-fld],
                    intercept = param$intercept,
                    est.b = param$csuv.mod.collection[fld.i, ]))
    }, param = list(x = x, y = y, intercept = intercept, flds = flds,
                    csuv.mod.collection = csuv.mod.collection))
    parallel::stopCluster(cl)
    rm(cl)
  } else{
    result <- lapply(seq_len(length(flds)),
                     function(fld.i) {
                       fld <- flds[[fld.i]]
                       return(lm.PR(x = x[-fld, ], y = y[-fld],
                                     intercept = intercept,
                                     est.b = csuv.mod.collection[fld.i, ]))
                     })
  }
  return(do.call(rbind, result))
}
