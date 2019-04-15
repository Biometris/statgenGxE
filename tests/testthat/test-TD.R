context("Class TD")

test_that("createTD creates objects of class TD", {
  expect_is(createTD(data = testData), "TD")
})

test_that("renaming columns 'one to one' works properly in createTD", {
  expect_equal(colnames(createTD(data = testData, genotype = "seed")[[1]])[1:2],
               c("genotype", "family"))
  expect_equal(colnames(createTD(data = testData, genotype = "seed",
                                 trial = "field")[[1]])[1:3],
               c("genotype", "family", "trial"))
  expect_error(createTD(data = testData, genotype = "a"), "has to be NULL or")
})

test_that("renaming columns 'one to many' works properly in createTD", {
  expect_equal(ncol(createTD(data = testData, rowId = "Y",
                             rowCoord = "Y")[[1]]), ncol(testData) + 1)
  expect_equal(colnames(createTD(data = testData, rowId = "Y",
                                 rowCoord = "Y")[[1]])[c(7, 13)],
  c("rowId", "rowCoord"))
})

test_that("class conversion works properly in createTD", {
  expect_is(createTD(data = testData)[[1]][, "Y"], "numeric")
  expect_is(createTD(data = testData, rowId = "Y")[[1]][, "rowId"], "factor")
  expect_is(createTD(data = testData,
                     rowCoord = "Y")[[1]][, "rowCoord"], "numeric")
  expect_is(createTD(data = testData)[[1]][, "checkId"], "factor")
})

test_that("attribute renamed is properly filled in createTD", {
  expect_null(attr(createTD(data = testData)[[1]], "renamed"))
  expect_is(attr(createTD(data = testData, genotype = "seed")[[1]], "renamed"),
            "data.frame")
  expect_equal(attr(createTD(data = testData, rowId = "Y",
                             rowCoord = "Y")[[1]], "renamed"),
               data.frame(orig = c("Y", "Y"),
                          new = c("rowId", "rowCoord"),
                          stringsAsFactors = FALSE))
})

test_that("dropTD functions properly", {
  testTD <- createTD(data = testData, trial = "field")
  expect_warning(dropTD(TD = testTD, rmTrials = "E4"),
                "The following trials are not in TD")
  expect_warning(dropTD(TD = testTD, rmTrials = c("E1", "E2", "E3")),
                 "All trials have been removed from TD")
  testTDrm <- dropTD(TD = testTD, rmTrials = "E1")
  expect_is(testTDrm, "TD")
  expect_named(testTDrm, c("E2", "E3"))
})

test_that("getMeta functions properly", {
  TD1 <- createTD(data = testData)
  meta1 <- getMeta(TD1)
  ## No trial defined, so only 1 row in meta
  expect_equal(nrow(meta1), 1)
  expect_equal(rownames(meta1), "testData")
  TD2 <- createTD(data = testData, trial = "field")
  meta2 <- getMeta(TD2)
  expect_equal(nrow(meta2), 3)
  expect_equal(rownames(meta2), c("E1", "E2", "E3"))
  expect_equal(meta2$trLocation, c("E1", "E2", "E3"))
})

test_that("setMeta functions properly", {
  TD1 <- createTD(data = testData, trial = "field")
  meta1 <- getMeta(TD1)
  meta1$trDesign <- c("res.rowcol", "rowcol", "res.ibd")
  TD2 <- setMeta(TD = TD1, meta = meta1)
  expect_equal(attr(x = TD2$E1, which = "trDesign"), "res.rowcol")
})

test_that("attribute design is properly filled in create TD", {
  expect_null(attr(createTD(data = testData)[[1]], "design"))
  expect_equal(attr(createTD(data = testData, trDesign = "rcbd")[[1]],
                    "trDesign"), "rcbd")
  expect_error(createTD(data = testData, trDesign = "abc"), "should be one of")
})

test_that("summary.TD produces correct output", {
  sumTD <- summary(createTD(data = testData), traits = c("t1", "t4"),
                   what = "all")
  expect_is(sumTD, "array")
  expect_equal(mean(testData$t1), sumTD["Mean", "t1", 1])
  expect_equal(max(testData$t4, na.rm = TRUE), sumTD["Max", "t4", 1])
})

test_that("option groupBy in summary.TD produces correct output", {
  sumTD <- summary(createTD(data = testData), traits = c("t1", "t4"),
                   groupBy = "field")
  expect_equal(dim(sumTD), c(3, 2, 3))
  expect_equivalent(as.numeric(by(data = testData$t1, INDICES = testData$field,
                                  FUN = mean)), sumTD["Mean", "t1", ])
  expect_equivalent(sumTD["Number of observations", "t4", ], c(27, 22, 26))
})

test_that("createTD accepts tibbles as input", {
  ## Skip on cran since it needs package tibble as extra dependency.
  skip_on_cran()
  expect_is(createTD(data = tibble::as.tibble(testData)), "TD")
})
