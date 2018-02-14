context("Class TD")

test_that("createTD creates objects of class TD", {
  expect_is(createTD(data = testData), "TD")
})

test_that("renaming columns 'one to one' works properly in createTD", {
  expect_equal(colnames(createTD(data = testData, genotype = "seed"))[1:2],
               c("genotype", "family"))
  expect_equal(colnames(createTD(data = testData, genotype = "seed",
                                 env = "field"))[1:3],
               c("genotype", "family", "env"))
  expect_error(createTD(data = testData, genotype = "a"), "has to be NULL or")
})

test_that("renaming columns 'one to many' works properly in createTD", {
  expect_equal(ncol(createTD(data = testData, rowId = "Y", rowCoordinates = "Y")),
               ncol(testData) + 1)
  expect_equal(colnames(createTD(data = testData, rowId = "Y",
                                 rowCoordinates = "Y"))[c(7, 13)],
  c("rowId", "rowCoordinates"))
})

test_that("class conversion works properly in createTD", {
  expect_is(createTD(data = testData)[, "Y"], "numeric")
  expect_is(createTD(data = testData, rowId = "Y")[, "rowId"], "factor")
  expect_is(createTD(data = testData,
                     rowCoordinates = "Y")[, "rowCoordinates"], "numeric")
  expect_is(createTD(data = testData)[, "checkId"], "factor")
})

test_that("attribute renamed is properly filled in createTD", {
  expect_null(attr(createTD(data = testData), "renamed"))
  expect_is(attr(createTD(data = testData, genotype = "seed"), "renamed"), "data.frame")
  expect_equal(attr(createTD(data = testData, rowId = "Y",
                             rowCoordinates = "Y"), "renamed"),
               data.frame(orig = c("Y", "Y"), new = c("rowId", "rowCoordinates"),
                          stringsAsFactors = FALSE))
})

test_that("attribute design is properly filled in create TD", {
  expect_null(attr(createTD(data = testData), "design"))
  expect_equal(attr(createTD(data = testData, design = "rcbd"), "design"), "rcbd")
  expect_error(createTD(data = testData, design = "abc"), "should be one of")
})

sumTD <- summary(createTD(data = testData), traits = c("t1", "t4"), what = "all")
test_that("summary.TD produces correct output", {
  expect_is(sumTD, "table")
  expect_equal(mean(testData$t1), sumTD["Mean", "t1"])
  expect_equal(max(testData$t4, na.rm = TRUE), sumTD["Max", "t4"])
})


