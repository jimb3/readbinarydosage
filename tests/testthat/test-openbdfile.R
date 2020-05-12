test_that("openbdfile", {
  bdfile <- system.file("extdata",
                        "bdinfo_set1a_1_1.rds",
                        package = "readbinarydosage")
  bdinfo <- readRDS(bdfile)
  snps <- readsnps(bdfile = bdinfo,
                   dosageonly = TRUE)
  expect_equal(snps$status, "Good")

  bdinfo$filename <- "xyz123abc456.bdose"
  snps <- readsnps(bdfile = bdinfo,
                   dosageonly = TRUE)
  expect_equal(snps$status, "Failed to open file")
  
  testfile <- system.file("extdata",
                          "error1.bdose",
                          package = "readbinarydosage")
  bdinfo$filename <- testfile
  snps <- readsnps(bdfile = bdinfo,
                   dosageonly = TRUE)
  expect_equal(snps$status, "Error reading binary dosage file header")
  
  testfile <- system.file("extdata",
                          "error2.bdose",
                          package = "readbinarydosage")
  bdinfo$filename <- testfile
  snps <- readsnps(bdfile = bdinfo,
                   dosageonly = TRUE)
  expect_equal(snps$status, "Not a binary dosage file")
  
  testfile <- system.file("extdata",
                          "error3.bdose",
                          package = "readbinarydosage")
  bdinfo$filename <- testfile
  snps <- readsnps(bdfile = bdinfo,
                   dosageonly = TRUE)
  expect_equal(snps$status, "Error reading format")
  
  testfile <- system.file("extdata",
                          "error4.bdose",
                          package = "readbinarydosage")
  bdinfo$filename <- testfile
  snps <- readsnps(bdfile = bdinfo,
                   dosageonly = TRUE)
  expect_equal(snps$status, "Unknown format")
  
  testfile <- system.file("extdata",
                          "error5.bdose",
                          package = "readbinarydosage")
  bdinfo$filename <- testfile
  snps <- readsnps(bdfile = bdinfo,
                   dosageonly = TRUE)
  expect_equal(snps$status, "Unknown subformat for format 1 or 2")
  
  testfile <- system.file("extdata",
                          "error6.bdose",
                          package = "readbinarydosage")
  bdinfo$filename <- testfile
  snps <- readsnps(bdfile = bdinfo,
                   dosageonly = TRUE)
  expect_equal(snps$status, "Unknown subformat for format 3 or 4")

  testfile <- system.file("extdata",
                          "test1.bdose",
                          package = "readbinarydosage")
  snps <- readsnps(bdfile = testfile,
                   dosageonly = TRUE)
  expect_equal(snps, "readsnps currently requires bdinfo as the input")
  
})
