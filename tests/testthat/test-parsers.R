suppressPackageStartupMessages(library(simrel2))
suppressPackageStartupMessages(library(testthat))

context("Testing all the parser functions.")


testthat::test_that(
  "Test char2list", {
    expect_equal(chr2list("1, 2, 3"), c(1, 2, 3))
    expect_equal(chr2list("1, 2; 3, 4"), list(c(1, 2), c(3, 4)))
    expect_equal(chr2list("1:2; 3:7"), list(1:2, 3:7))
    expect_equal(chr2list("1, 2, 3;"), list(c(1, 2, 3)))
  }
)

testthat::test_that(
  "Test list2chr", {
    expect_equal(list2chr(c(1, 2, 3)), "1; 2; 3")
    expect_equal(list2chr(list(c(1, 2), c(3, 4))), "1, 2; 3, 4")
    expect_equal(list2chr(list(1:2, 3:7)), "1, 2; 3, 4, 5, 6, 7")
  }
)
