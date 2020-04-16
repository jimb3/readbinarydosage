#' @useDynLib readbinarydosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

release_questions <- function() {
  c("Are you sure you want to do this?",
    "Are you really sure?")
}

#' readsnps
#' 
#' Test is the file can be opened
#' 
#' @param bdfile Name of binary dosage data file
#' 
#' @return
#' Returns status of reading file
#' @export
#'
#' @examples
#' testfile <- system.file("extdata",
#'                         "test1.bdose",
#'                         package = "readbinarydosage")
#' readsnps(bdfile = testfile)
readsnps <- function(bdfile) {
  return(readsnpsc(filenames = bdfile))
}