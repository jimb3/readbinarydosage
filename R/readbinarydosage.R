#' @useDynLib readbinarydosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

release_questions <- function() {
  c("Are you sure you want to do this?",
    "Are you really sure?")
}

# Create blocks used for faster reading of binary dosage files
assignblocks <- function(nsub, nsnps, indices, dsize) {
  if (nsub < 10000) {
    blksize <- 5000
  } else if (nsub < 25000) {
    blksize <- 2000
  } else if (nsub < 50000) {
    blksize <- 1000
  } else if (nsub < 100000) {
    blksize <- 500
  } else {
    blksize <- 200
  }
  
  nblks <- ceiling(nsnps / blksize)
  if (nblks == 1) {
    blksize <- nsnps
    fsnp <- 1
    blkbytes <- indices[nsnps] - indices[1] + dsize[nsnps]
  } else {
    fsnp <- seq(1, (nblks - 1) * blksize + 1, blksize)
    blkbytes <- numeric(nblks)
    blkbytes[1:(nblks - 1)] <- indices[fsnp[2:nblks]] - indices[fsnp[1:(nblks - 1)]]
    blkbytes[nblks] <- indices[nsnps] - indices[fsnp[nblks]] + dsize[nsnps]
  }
  return(list(
    blksize = blksize,
    nblks = nblks,
    fsnp = fsnp,
    blkbytes = blkbytes))
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
  if (is.character(bdfile)) {
    print("readsnps currently requires bdinfo as the input")
    return ("readsnps currently requires bdinfo as the input")
  }
  if (class(bdfile) == FALSE) {
    print("bdfile not of type character or bdose-info")
    return ("bdfile not of type character or bdose-info")
  }
  blocks <- assignblocks(nsub = nrow(bdfile$samples),
                         nsnps = length(bdfile$snps$snpid),
                         indices = bdfile$indices,
                         dsize = bdfile$datasize)
  return(readsnpsc(filename = bdfile$filename,
                   indices = bdfile$indices,
                   numsub = nrow(bdfile$samples),
                   TRUE))
}