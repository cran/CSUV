get.fit.methods <- function() {
  fit.methods <- c(lasso = lm.lasso.min,
                   elastic = lm.elastic.half,
                   relaxo = lm.relaxo,
                   mcp = lm.mcp,
                   scad = lm.scad)
  return(fit.methods.helper.add.name(fit.methods))
}

get.path.fit.methods <- function() {
  path.fit.methods <- list(lasso = lm.lasso.path,
                           elastic = lm.elastic.half.path,
                           relaxo = lm.relaxo.path,
                           scad = lm.scad.path,
                           mcp = lm.mcp.path)
  return(fit.methods.helper.add.name(path.fit.methods))
}

## ----- helper function -----
fit.methods.helper.add.name <- function(fit.methods) {
  for (fit.method.name in names(fit.methods)) {
    fit.method <- fit.methods[[fit.method.name]]
    attr(fit.method, "name") <- fit.method.name
    fit.methods[[fit.method.name]] <- fit.method
  }
  return(fit.methods)
}
