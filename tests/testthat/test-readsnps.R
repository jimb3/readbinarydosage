test_that("readsnps", {
  bdfile <- system.file("extdata",
                        "test1.bdose",
                        package = "readbinarydosage")
  expect_equal(readsnps(bdfile = bdfile, TRUE),
               "readsnps currently requires bdinfo as the input")

  bdfile <- system.file("extdata",
                        "bdinfo_set1a_1_1.rds",
                        package = "readbinarydosage")
  bdinfo <- readRDS(bdfile)
  
  expect_equal(readsnps(bdfile = bdinfo, 1L),
               "dosageonly must be a logical value")

  snps <- readsnps(bdfile = bdinfo, FALSE)
  expect_equal(snps$status, "Good")
})
