test_that("readsnps", {
  bdfile <- system.file("extdata",
                        "test1.bdose",
                        package = "readbinarydosage")
  snps <- readsnps(bdfile = bdfile)
  expect_equal(snps, "readsnps currently requires bdinfo as the input")

  bdfile <- system.file("extdata",
                        "bdinfo_set1a_1_1.rds",
                        package = "readbinarydosage")
  bdinfo <- readRDS(bdfile)
  snps <- readsnps(bdfile = bdinfo)
  expect_equal(snps$status, "Good")
})
