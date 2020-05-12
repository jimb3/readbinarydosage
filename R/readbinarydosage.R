#' @useDynLib readbinarydosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

release_questions <- function() {
  c("Are you sure you want to do this?",
    "Are you really sure?")
}

# Create blocks used for faster reading of binary dosage files
assignblocks <- function(nsub, nsnps, snploc, snpbytes) {
  if (nsub < 10000) {
    blksnps <- 5000
  } else if (nsub < 25000) {
    blksnps <- 2000
  } else if (nsub < 50000) {
    blksnps <- 1000
  } else if (nsub < 100000) {
    blksnps <- 500
  } else if (nsub < 250000) {
    blksnps <- 200
  } else if (nsub < 500000) {
    blksnps <- 100
  } else {
    blksnps <- 50
  }
  
  nblks <- ceiling(nsnps / blksnps)
  if (nblks == 1) {
    blksnps <- nsnps
    fsnp <- 1
    blkloc <- snploc[fsnp]
    blkbytes <- snploc[nsnps] - snploc[1] + snpbytes[length(snpbytes)]
  } else {
    fsnp <- seq(1, (nblks - 1) * blksnps + 1, blksnps)
    blkloc <- snploc[fsnp]
    blkbytes <- numeric(nblks)
    blkbytes[1:(nblks - 1)] <- snploc[fsnp[2:nblks]] - snploc[fsnp[1:(nblks - 1)]]
    blkbytes[nblks] <- snploc[nsnps] - snploc[fsnp[nblks]] + snpbytes[length(snpbytes)]
    blksnps <- rep(blksnps, nblks)
    blksnps[nblks] <- nsnps %% blksnps[1]
  }
  return(list(
    nblks = nblks,
    blksnps = blksnps,
    fsnp = fsnp,
    blkloc = blkloc,
    blkbytes = blkbytes))
}

#' readsnps
#' 
#' Test is the file can be opened
#' 
#' @param bdfile Name of binary dosage data file
#' @param subjects List of indices or subjects to get values for
#' @param dosageonly indicator if dosage values only are to be returned
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
readsnps <- function(bdfile,
                     subjects,
                     snps,
                     dosageonly) {
  if (is.logical(dosageonly) == FALSE) {
    print("dosageonly must be a logical value")
    return("dosageonly must be a logical value")
  }
  if (is.character(bdfile)) {
    print("readsnps currently requires bdinfo as the input")
    return ("readsnps currently requires bdinfo as the input")
  }
  if (class(bdfile) == FALSE) {
    print("bdfile not of type character or bdose-info")
    return ("bdfile not of type character or bdose-info")
  }
  if (missing(subjects) == TRUE) {
    subjects <- 1:nrow(bdfile$samples)
  } else {
    if (is.numeric(subjects) == FALSE)
      return ("subjects must be a numeric vector")
    if (all(as.integer(subjects) == subjects) == FALSE)
      return ("subjects must be a integer vector")
  }
  if (missing(snps) == TRUE) {
    snps <- 1:length(bdfile$snps$snpid)
  } else {
    if (is.numeric(snps) == FALSE)
      return ("snps must be a numeric vector")
    if (all(as.integer(snps) == snps) == FALSE)
      return ("snps must be a integer vector")
  }
  blocks <- assignblocks(nsub = nrow(bdfile$samples),
                         nsnps = length(bdfile$snps$snpid),
                         snploc = bdfile$indices,
                         snpbytes = bdfile$datasize)
#  return(blocks)
  return(readsnpsc(filename = bdfile$filename,
                   subjects = subjects,
                   snps = snps,
                   nsub = nrow(bdfile$samples),
                   snploc = bdfile$indices,
                   snpbytes = bdfile$datasize,
                   blksnps = blocks$blksnps,
                   firstsnp = blocks$fsnp,
                   blkloc = blocks$blkloc,
                   blkbytes = blocks$blkbytes, 
                   dosageonly))
}