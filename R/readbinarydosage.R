#' @useDynLib readbinarydosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' readsnps
#' Does nothing
#' @return
#' Returns 1
#' @export
#'
#' @examples
#' readsnps()
readsnps <- function() {
  return(readsnpsc())
}