###########################################################
## This file includes functions related to plotting for CSUV
## (include helper functions)
###########################################################

# @import graphics

## ------- plot helper functions ---------
#' Graphical illustration of selection uncertainty
#' @export
#' @param x fitted results from CSUV::csuv()
#' @param with.unconditional TRUE to get a unconditonal boxplot on the same graph. Default is FALSE
#' @param compare.method.fit (optional) fitted results from CSUV::lm.compare.methods(). Alternatively, user can provide a data frame with each row contains the estimated coefficients from a method. The name of each row should be corresponding to the name of the method. The first value of each row should be the value of the intercept
#' @param cv.mod (optional) a vector of estimated coefficients from cross validation. The first value should be the value of the intercept
#' @param with.thr whether the selection by the CSUV should be show. Default is TRUE
#' @param with.violin whether the graph with violin plot
#' @param to.shade whether to shade the graph by the relative frequency calculated by CSUV. Default is TRUE
#' @param ci.method how the confidence interval should be calculated. Default is "conditional"
#' @param level the significant level of the whiskers. Default is 0.1
#' @param var.freq.thr minimum variable frequency to show, default is 0.1
#' @param log.level log level to set. Default is NULL, which means no change in log level. See the function CSUV::set.log.level for more details
#' @param ... additional argument for plot
#' @return a ggplot object
#' @examples
#' \donttest{
#' X = matrix(rnorm(1000), nrow = 100)
#' Y = rowSums(X[,1:3])+rnorm(100)
#' mod.0 = csuv(X, Y, intercept = FALSE, q = 0, method.names = NULL)
#' cv.mod = lm.cv(X, Y, intercept = FALSE, fit.percent = 0.5, num.repeat = 50)
#' compare.mod = lm.compare.method(X, Y, intercept = FALSE)
#' plot(mod.0, compare.method.fit = compare.mod, cv.mod = cv.mod$est.b)
#' }
plot.csuv <- function(x,
                    with.unconditional = FALSE,
                    compare.method.fit = NULL,
                    cv.mod = NULL,
                    with.thr = TRUE,
                    with.violin = FALSE,
                    to.shade = TRUE,
                    ci.method = "conditional",
                    level = 0.1,
                    var.freq.thr = 0.1,
                    log.level = NULL,
                    ...) {
  # plot.new()
  graphics::par(new = FALSE)
  if (!is.null(log.level)){
    set.log.level(log.level)
  }
  return(csuv.plot.helper(new.fit = x,
                           with.unconditional = with.unconditional,
                           compare.method.fit = compare.method.fit,
                           compare.method.names = rownames(compare.method.fit),
                           cv.mod = cv.mod,
                           ci.method = ci.method,
                           print.compare.method.points = FALSE,
                           with.thr = with.thr,
                           with.violin = with.violin,
                           to.shade = to.shade,
                           level = level,
                           var.freq.thr = var.freq.thr, ...))
}


## ------- plot helper functions ---------
#' Helper function, please do not use it
#' @export csuv.plot.helper
#' @param new.fit fitted results from CSUV::csuv()
#' @param with.unconditional TRUE to get a unconditonal boxplot on the same graph. Default is FALSE
#' @param compare.method.fit (optional) fitted results from CSUV::lm.compare.methods()
#' @param compare.method.names (optional) names of method to compare
#' @param cv.mod (optional) fitted results from cross validation
#' @param print.compare.method.points Default is FALSE
#' @param with.thr whether the selection by the CSUV should be show. Default is TRUE
#' @param with.violin whether the graph with violin plot
#' @param to.shade whether to shade the graph by the relative frequency calculated by CSUV. Default is TRUE
#' @param ci.method how the confidence interval should be calculated. Default is "conditional"
#' @param level the significant level of the whiskers. Default is 0.1
#' @param var.freq.thr minimum variable frequency to show, default is 0.1
#' @param ... additional argument for plot
#' @return a ggplot object
csuv.plot.helper <- function(new.fit,
                             with.unconditional = FALSE,
                             compare.method.fit = NULL,
                             compare.method.names = NULL,
                             cv.mod = NULL,
                             print.compare.method.points = FALSE,
                             ci.method = "conditional",
                             with.thr = TRUE,
                             with.violin = FALSE,
                             to.shade = TRUE,
                             level = 0.1,
                             var.freq.thr = 0.1,
                             ...) {
  shiny::req(new.fit)

  ## Just to quiet CRAN
  variables <- NULL
  type <- NULL
  coefficients <-NULL
  cutoff <- NULL

  mod <- new.fit$mod.collection
  if (is.null(colnames(mod))) { # make sure mod has colnames
    colnames(mod) <- paste0("X", 0:(ncol(mod) - 1))
  }
  original.method.names <- names(get.compare.methods())
  if (length(compare.method.fit) && is.null(compare.method.names)) {
    compare.method.names <- rownames(compare.method.fit)
  }
  compare.method.fit <- compare.method.fit[compare.method.names, , drop = FALSE]


  ## ======== get the var.order  =======
  var.freq <- new.fit$variable.freq
  var.order <- get.plot.var.order(new.fit$variable.order, var.freq,
                                  compare.method.fit, cv.mod, var.freq.thr = var.freq.thr)

  ## ======== convert to ggplot df  =======
  ggplot.df <- get.df.for.gg.plot(mod = mod[, -1], csuv.mod = new.fit$est.b["csuv.m", -1],
                                  compare.mod = compare.method.fit[, -1, drop = FALSE],
                                  cv.mod = cv.mod[-1],
                                  tau = var.freq,
                                  var.order = var.order,
                                  print.compare.method.points = print.compare.method.points)

  ## ====== prepare info for plotting =======
  shade.info <- get.shade.col.n.lab(var.freq)
  get.cond.boxplot.stat <- function(x) {
    get.boxplot.stat(x, is.conditional = TRUE)
  }
  get.cond.errorbar.stat <- function(x) {
    get.errorbar.stat(x, is.conditional = TRUE)
  }

  get.boxplot.stat <- function(x, is.conditional = FALSE) {
    r <- rep(0, 6)
    cond.x <- x
    if (is.conditional) cond.x <- x[which(x != 0)]
    if (length(cond.x)) {
      tau <- ifelse(is.conditional, max(mean(x > 0), mean(x < 0)) * 0.8, 0.8)
      r <- c(stats::quantile(cond.x, probs = c(0.25, 0.25, 0.5, 0.75, 0.75)), tau)
    }
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax", "width")
    return(r)
  }
  get.errorbar.stat <- function(x, is.conditional = FALSE) {
    r <- rep(0, 3)
    cond.x <- x
    if (is.conditional) cond.x <- x[which(x != 0)]
    if (length(cond.x)) {
      tau <- ifelse(is.conditional, max(mean(x > 0), mean(x < 0)) * 0.8, 0.8)
      r <- c(stats::quantile(cond.x, probs = c((level / 2), (1 - level / 2))), tau)
    }
    names(r) <- c("ymin", "ymax", "width")
    return(r)
  }

  y.min <- min(ggplot.df[, "coefficients"])
  method.point.col <- c(csuv_m = "red")
  method.point.sty <- c(csuv_m = 16)
  if (!is.null(cv.mod)) {
    method.point.col <- c(method.point.col, cv = "blue")
    method.point.sty <- c(method.point.sty, cv = 1)
  }
  if (print.compare.method.points & !is.null(compare.method.fit)) {
    col <- grDevices::rainbow(length(original.method.names) + 1)[-1]
    sty <- rep(16, length(original.method.names))
    names(col) <- original.method.names
    names(sty) <- original.method.names
    method.point.col <- c(method.point.col, col)
    method.point.sty <- c(method.point.sty, sty)
  }
  # if (with.thr) {
  #   method.point.col = c(method.point.col, csuv_m = "green", csuv_s = "yellow")
  # }
  ## ========= plot ============
  p <- ggplot2::ggplot(ggplot.df, ggplot2::aes(x = variables, y = coefficients))

  if (to.shade) {
    p <- p + ggplot2::geom_tile(data = subset(ggplot.df, type == "tau"),
                                ggplot2::aes(fill = factor(pmin(coefficients %/% 10 * 10, 90),
                                                           levels = shade.info$tau.factors, ordered = TRUE),
                                             y = 0, height = Inf),
                                alpha = 0.3) +
      ggplot2::scale_fill_manual(name = bquote(tau), values = shade.info$col.rgb,
                        label = shade.info$col.labels)
  }
  p <- p + ggplot2::geom_hline(yintercept = 0)


  # == add the box plot ==
  p <- p + ggplot2::stat_summary(data = subset(ggplot.df, type == "fit"),
                                fun.data = get.cond.errorbar.stat, geom = "errorbar", col = "black", lwd = 0.7) +
    ggplot2::stat_summary(data = subset(ggplot.df, type == "fit"),
                          fun.data = get.cond.boxplot.stat, geom = "boxplot", varwidth = TRUE, fill = "gold", color = "orangered")
  if (with.unconditional) {
    p <- p + ggplot2::stat_summary(data = subset(ggplot.df, type == "fit"),
                                   fun.data = get.errorbar.stat, geom = "errorbar", col = "brown", lwd = 0.7) +
      ggplot2::stat_summary(data = subset(ggplot.df, type == "fit"),
                            fun.data = get.boxplot.stat, geom = "boxplot", varwidth = TRUE, fill = "palegreen", color = "royalblue", alpha = 0.5)
  }

  # == add the violin plot ==
  if (with.violin) {
    p <- p + ggplot2::geom_violin(data = subset(ggplot.df, type == "fit" & coefficients != 0), alpha = 0.5, col = "blue", scale = "count")
  }

  # == add coef estimation points ==
  p <- p + ggplot2::geom_point(data = subset(ggplot.df, type %in% names(method.point.col)),
                               ggplot2::aes(color = type, shape = type)) +
    ggplot2::scale_shape_manual(name = "methods", values = method.point.sty) +
    ggplot2::scale_color_manual(name = "methods", values = method.point.col)

  # == add tau and selection frequency text ==
  p <- p + ggplot2::geom_text(data = subset(ggplot.df, type == "tau"),
                              ggplot2::aes(y = y.min - 1, label = coefficients)) +
    ggplot2::geom_text(label = "tau", x = 1, y = y.min - 0.75, parse = TRUE,
                       hjust = 0.5, vjust = 0.5) #+
    # ggplot2::coord_cartesian(xlim = c(1, x.max+0.5), clip = "off")

  if (!is.null(compare.method.fit)) {
    p <- p + ggplot2::geom_text(data = subset(ggplot.df, type == "compare"),
                                ggplot2::aes(y = y.min - 1.5, label = coefficients), col = "blue") +
      ggplot2::geom_text(label = "comparing methods", x = 0.6, y = y.min - 1.25,
                         col = "blue", hjust = 0, vjust = 0.5)
  }

  # == add threshold ==
  cf <- data.frame("type" = c("csuv_m", "csuv_s"),
                   cutoff = c((sum(new.fit$est.b["csuv.m", ] != 0) + 0.5),
                              (sum(new.fit$est.b["csuv.s", ] != 0) + 0.5)), dummy = 0)
  if (with.thr) {
    p <- p + ggplot2::geom_vline(data = cf,
                                 ggplot2::aes(xintercept = cutoff, linetype = type),
                                 color = c("chartreuse2", "blue"), size = c(1.5, 1)) +
      ggplot2::scale_linetype_manual(name = "cutoff",
                                     values = c(csuv_m = "solid", csuv_s = "dashed"),
                                     guide = ggplot2::guide_legend(override.aes = list(color = c("chartreuse2", "blue"), size = c(1.5, 1))))
  }

  # update axis
  p <- p + ggplot2::ylab("estimated coefficients") +
    ggplot2::xlab("sorted covariates") +
    ggplot2::theme_bw()
  return(p)
}
