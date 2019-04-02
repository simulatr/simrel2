suppressPackageStartupMessages(library(simrel2))
suppressPackageStartupMessages(library(testthat))

context("Testing all the utilities functions.")


testthat::test_that(
  "Test depth of a list", {
    expect_equal(depth(list(list(1, 2), list(2, 3))), 3)
    expect_equal(depth(1:5), 1)
    expect_equal(depth(list(list(1, 2), list(2, list(3, 4)))), 4)
    expect_equal(depth(list(1:2, 3:4)), 2)
  }
)

