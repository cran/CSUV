#' Set the level of logger
#' @export set.log.level
#' @param level log level, setting the level to futile.logger::DEBUG provides most details log, whereas setting the level to futile.logger::WARN provides least details log
#' @return None
#' @examples
#' \donttest{
#' set.log.level(futile.logger::DEBUG)
#' set.log.level(futile.logger::INFO)
#' set.log.level(futile.logger::WARN)
#' }
set.log.level <- function(level) {
  futile.logger::flog.threshold(level)
}
