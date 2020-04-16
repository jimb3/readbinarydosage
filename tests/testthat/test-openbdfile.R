test_that("openbdfile", {
  errorfile <- "xyz123abc456.bdose"
  snps <- readsnps(bdfile = errorfile)
  expect_equal(snps$status, "Failed to open file")
  
  errorfile <- system.file("extdata",
                           "error1.bdose",
                           package = "readbinarydosage")
  snps <- readsnps(bdfile = errorfile)
  expect_equal(snps$status, "Error reading binary dosage file header")
  
  errorfile <- system.file("extdata",
                           "error2.bdose",
                           package = "readbinarydosage")
  snps <- readsnps(bdfile = errorfile)
  expect_equal(snps$status, "Not a binary dosage file")
  
  errorfile <- system.file("extdata",
                           "error3.bdose",
                           package = "readbinarydosage")
  snps <- readsnps(bdfile = errorfile)
  expect_equal(snps$status, "Error reading format")
  
  errorfile <- system.file("extdata",
                           "error4.bdose",
                           package = "readbinarydosage")
  snps <- readsnps(bdfile = errorfile)
  expect_equal(snps$status, "Unknown format")
  
  errorfile <- system.file("extdata",
                           "error5.bdose",
                           package = "readbinarydosage")
  snps <- readsnps(bdfile = errorfile)
  expect_equal(snps$status, "Unknown subformat for format 1 or 2")
  
  errorfile <- system.file("extdata",
                           "error6.bdose",
                           package = "readbinarydosage")
  snps <- readsnps(bdfile = errorfile)
  expect_equal(snps$status, "Unknown subformat for format 3 or 4")

  testfile <- system.file("extdata",
                          "test1.bdose",
                          package = "readbinarydosage")
  snps <- readsnps(bdfile = testfile)
  expect_equal(snps$status, "Good")
})
