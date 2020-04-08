test_that("readsnps", {
  snps <- readsnps()
  expect_equal(snps$x, 1)
})
