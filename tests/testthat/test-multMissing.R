context("multMissing")

D <- data.frame(c1 = 1:5, c2 = c(NA, 1, 2, NA, 3))
M1 <- matrix(c(1, 2, 3, NA, NA, 5, 6,
               6, 5, 1, NA, 3, 2, 1), nrow = 7, ncol = 2)
M2 <- matrix(c("1", "2", "3", NA, "b", "5", "6",
               "6", "5", "b", NA, "3", "2", "1"), nrow = 7, ncol = 2)
v <- c(NA, 1, 2, NA, 3, 4, 5)

test_that("general checks in multMissing function properly", {
  expect_error(multMissing(NULL),
               "Y should be a data.frame, matrix or vector")
  expect_error(multMissing(D, maxIter = 0),
               "a single numerical value greater than or equal to 1")
  expect_error(multMissing(D, naStrings = 1),
               "NULL or a character vector")
})

test_that("multMissing preserves input class", {
  expect_is(multMissing(D), "data.frame")
  expect_is(multMissing(M1), "matrix")
  expect_is(multMissing(v), "numeric")
})

test_that("multMissing leaves input without missings intact", {
  expect_warning(v1 <- multMissing(1:3), "no missing values found")
  expect_equal(v1, 1:3)
  expect_warning(m1 <- multMissing(matrix(1:4, nrow = 2)),
                 "no missing values found")
  expect_equivalent(m1, matrix(1:4, nrow = 2))
})

test_that("multMissing imputes correct values", {
  expect_equal(multMissing(D)[["c2"]], c(0.5, 1, 2, 2.42857142857143, 3))
  expect_equal(multMissing(M1)[, 1], c(1, 2, 3, 3.4, 3.4, 5, 6))
  expect_equal(multMissing(M1)[, 2], c(6, 5, 1, 3, 3, 2, 1))
  expect_equal(multMissing(v), c(3, 1, 2, 3, 3, 4, 5))
})

test_that("option maxIter functions properly", {
  expect_warning(mmD1 <- multMissing(D, maxIter = 1), "No convergence achieved")
  expect_equal(mmD1[["c2"]], c(0.5, 1, 2, 2.42857142857143, 3))
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

test_that("matrix with row with only NA cannot be imputed", {
          expect_error(multMissing(matrix(c(NA, NA, 1:2), nrow = 2)),
                       "At least one unit must have no missing values")
})
