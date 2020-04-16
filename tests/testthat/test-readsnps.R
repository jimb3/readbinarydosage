test_that("readsnps", {
  testfile <- system.file("extdata",
                          "test1.bdose",
                          package = "readbinarydosage")
  snps <- readsnps(bdfile = testfile)
  expect_equal(snps$status, "Good")
})
