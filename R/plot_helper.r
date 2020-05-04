get.plot.var.order <- function(var.order, var.freq, compare.method.fit,
                               cv.fit, var.freq.thr) {
  if (is.null(compare.method.fit) & is.null(cv.fit)) {
    return(var.order[1:sum(var.freq >= var.freq.thr)])
  } else{
    cv.sel.i <- NULL
    compare.sel.i <- NULL
    if (!is.null(compare.method.fit)) compare.sel.i <- which(colSums(compare.method.fit[, -1, drop = FALSE] != 0) != 0)
    if (!is.null(cv.fit)) cv.sel.i <- which(var.freq < var.freq.thr & cv.fit[-1] != 0)
    low.tau.but.compare.i <- intersect(which(var.freq < var.freq.thr), union(compare.sel.i, cv.sel.i))

    if (length(low.tau.but.compare.i)) {
      cv.dummy = compare.dummy = rep(0, length(low.tau.but.compare.i))
      if (!is.null(cv.fit)) cv.dummy <- cv.fit[-1][low.tau.but.compare.i]
      if (!is.null(compare.method.fit)) {
        compare.dummy <- colMeans(compare.method.fit[, -1, drop = FALSE][, low.tau.but.compare.i, drop = FALSE] != 0)
      }
      var.order <- c(var.order[1:sum(var.freq >= var.freq.thr)],
                     low.tau.but.compare.i[order(-var.freq[low.tau.but.compare.i],
                                                 -abs(compare.dummy), -abs(cv.dummy))])
    }
  }
  return(var.order)
}

get.df.for.gg.plot <- function(mod, csuv.mod, compare.mod, cv.mod, tau,
                               var.order, print.compare.method.points) {
  # the mod, compare.mod, cv.mod all should have NO intercept
  long.df <- reshape2::melt(data = mod,
                            variables.name = "variables",
                            value.name = "coefficients"
  )
  long.df <- long.df[, -1] # remove the first col
  colnames(long.df)[1] <- "variables"
  long.df[, "type"] <- "fit"

  # ===== add other values ==========
  var.names <- colnames(mod)
  long.df <- rbind(long.df, data.frame(variables = var.names,
                                       coefficients = csuv.mod,
                                       type = "csuv_m"))
  long.df <- rbind(long.df, data.frame(variables = var.names,
                                       coefficients = round(tau * 100),
                                       type = "tau"))
  if (!is.null(compare.mod)) {
    long.df <- rbind(long.df, data.frame(variables = var.names,
                                         coefficients = round(colMeans(compare.mod != 0) * 100), type = "compare"))
    if (print.compare.method.points) {
      long.compare.df <- reshape2::melt(data = compare.mod,
                                        variables.name = "variables",
                                        value.name = "coefficients")
      colnames(long.compare.df)[1:2] <- c("type", "variables")
      long.df <- rbind(long.df, long.compare.df[, c("variables", "coefficients", "type")])
    }
  }
  if (!is.null(cv.mod)) long.df <- rbind(long.df, data.frame(variables = var.names,
                                                             coefficients = cv.mod,
                                                             type = "cv"))

  ordered.selected.variables <- var.names[var.order]
  long.df <- long.df[which(long.df[, "variables"] %in% ordered.selected.variables), ]
  long.df$variables <- factor(long.df$variables,
                              levels = ordered.selected.variables, ordered = TRUE)
  return(long.df)
}

get.shade.col.n.lab <- function(tau) {
  # tau is \in [0,1]
  tau.i <- as.numeric(sort(unique(pmin((tau * 100) %/% 10, 9)),
                           decreasing = TRUE))
  tau.factors <- as.factor(c((9:0) * 10))
  col.rgb <- sapply(tau.i, function(i) grDevices::rgb(1 - i / 10, 1 - i / 10, 1 - i / 10))
  col.labels <- sapply(tau.i, function(i)
    ifelse(i == 9, paste0("[", i * 10, ",", (i + 1) * 10, "]"),
           paste0("[", i * 10, ",", (i + 1) * 10, ")")))
  return(list(tau.factors = tau.factors,
               col.rgb = col.rgb,
               col.labels = col.labels))
}
