context("Report")

## AMMI

test_that("Test that createReport function only accepts correctly named output", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  geAmmi <- gxeAmmi(TD = BLUEs, trait = "t1")
  expect_error(report(geAmmi, outfile = tempfile(fileext = ".pd")),
               "Invalid output filename")
  expect_error(report(geAmmi, outfile = tempfile(fileext = "a b.pdf")),
               "outfile path cannot contain spaces")
  expect_silent(report(geAmmi, outfile = tempfile(fileext = ".pdf")))
})

## Finlay Wilkinson

test_that("FW report functions correctly", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  geFw <- gxeFw(TD = BLUEs, trait = "t1")
  expect_silent(report(geFw))
})

## Stability

test_that("Stability report functions correctly", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  geStab <- gxeStability(TD = BLUEs, trait = "t1")
  expect_silent(report(geStab))
})

## varCov

test_that("varCov report functions correctly", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  geVarCov <- gxeVarCov(TD = BLUEs, trait = "t1")
  expect_silent(report(geVarCov))
})

