context("Statistical Helper Functions")

test_that("function seVar gives correct results", {
  expect_equal(seVar(testData$t1), 102.714733134436)
  expect_equal(unname(seVar(testData[c("t1", "t2")])),
               c(102.714733134436, 0.0301924400875601))
  expect_equal(unname(seVar(as.matrix(testData[c("t1", "t2")]))),
               c(102.714733134436, 0.0301924400875601))
  expect_equal(seVar(testData$t3), NA_integer_)
  expect_equal(seVar(testData$t3, na.rm = TRUE), 16.5770418740614)
  expect_equal(unname(seVar(testData[c("t1", "t3")])),
               c(102.714733134436, NA))
  expect_equal(unname(seVar(testData[c("t1", "t3")], na.rm = TRUE)),
               c(102.714733134436, 16.5770418740614))
})

test_that("function skewness gives correct results", {
  expect_equal(skewness(testData$t1), 0.570922732616638)
  expect_equal(unname(skewness(testData[c("t1", "t2")])),
               c(0.570922732616638, 0.619398806071279))
  expect_equal(unname(skewness(as.matrix(testData[c("t1", "t2")]))),
               c(0.570922732616638, 0.619398806071279))
  expect_equal(skewness(testData$t3), NA_integer_)
  expect_equal(skewness(testData$t3, na.rm = TRUE), -0.221167185848707)
  expect_equal(unname(skewness(testData[c("t1", "t3")])),
               c(0.570922732616638, NA))
  expect_equal(unname(skewness(testData[c("t1", "t3")], na.rm = TRUE)),
               c(0.570922732616638, -0.221167185848707))
})

test_that("function seKurtosis gives correct results", {
  expect_equal(seSkewness(100), 0.24137977904013)
  expect_warning(seSkew <- seSkewness(2),
                 "the standard error of skewness cannot be calculated")
  expect_equal(seSkew, NA)
})

test_that("function kurtosis gives correct results", {
  expect_equal(kurtosis(testData$t1), 0.43414980172502)
  expect_equal(unname(kurtosis(testData[c("t1", "t2")])),
               c(0.434149801725015, -0.357971292973335))
  expect_equal(unname(kurtosis(as.matrix(testData[c("t1", "t2")]))),
               c(0.434149801725015, -0.357971292973335))
  expect_equal(kurtosis(testData$t3), NA_integer_)
  expect_equal(kurtosis(testData$t3, na.rm = TRUE), -0.340117318912085)
  expect_equal(unname(kurtosis(testData[c("t1", "t3")])),
               c(0.434149801725015, NA))
  expect_equal(unname(kurtosis(testData[c("t1", "t3")], na.rm = TRUE)),
               c(0.434149801725015, -0.340117318912085))
})

test_that("function seKurtosis gives correct results", {
  expect_equal(seKurtosis(100), 0.478331132994813)
  expect_warning(seKurt <- seKurtosis(2),
                 "the standard error of kurtosis cannot be calculated")
  expect_equal(seKurt, NA)
})


