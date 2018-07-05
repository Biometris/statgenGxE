context("multMissing")

D <- data.frame(c1 = 1:5, c2 = c(NA, 1, 2, NA, 3))
M1 <- matrix(c(1, 2, 3, NA, NA, 5, 6,
               6, 5, 1, NA, 3, 2, 1), nrow = 7, ncol = 2)
M2 <- matrix(c("1", "2", "3", NA, "b", "5", "6",
               "6", "5", "b", NA, "3", "2", "1"), nrow = 7, ncol = 2)
v <- c(NA, 1, 2, NA, 3, 4, 5)

test_that("multMissing preserves input class", {
  expect_is(multMissing(D), "data.frame")
  expect_is(multMissing(M1), "matrix")
  expect_is(multMissing(v), "numeric")
})

test_that("multMissing imputes correct values", {
  expect_equal(multMissing(D)$c2, c(0.5, 1, 2, 2.42857142857143, 3))
  expect_equal(multMissing(M1)[, 1], c(1, 2, 3, 3.4, 3.4, 5, 6))
  expect_equal(multMissing(M1)[, 2], c(6, 5, 1, 3, 3, 2, 1))
  expect_equal(multMissing(v), c(3, 1, 2, 3, 3, 4, 5))
})

test_that("option maxIter functions properly", {
  expect_warning(mmD1 <- multMissing(D, maxIter = 1), "No convergence achieved")
  expect_equal(mmD1$c2, c(0.5, 1, 2, 2.42857142857143, 3))
})

test_that("option naStrings functions properly", {
  ## No effect when strings don't occur.
  expect_equal(multMissing(M1, naStrings = c("a", "b"))[, 1],
               c(1, 2, 3, 3.4, 3.4, 5, 6))
  ## Warning when strings are not converted.
  expect_warning(multMissing(M2), "Warning when converting data to numeric")
  ## "b" treated as NA.
  expect_equal(multMissing(M2, naStrings = "b")[, 1],
               c(1, 2, 3, 3.49997515882245, 3.99985095293469, 5, 6))
  expect_equal(multMissing(M2, naStrings = "b")[, 2],
               c(6, 5, 3.99997226960459, 3.49999537826743, 3, 2, 1))
})
