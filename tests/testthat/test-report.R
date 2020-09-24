context("Report")

## AMMI

test_that("Test that createReport function only accepts correctly named output", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  skip_on_ci()
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
  skip_on_ci()
  geFw <- gxeFw(TD = BLUEs, trait = "t1", maxIter = 30)
  expect_silent(report(geFw, outfile = tempfile(fileext = ".pdf")))
})

## Stability

test_that("Stability report functions correctly", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  skip_on_ci()
  geStab <- gxeStability(TD = BLUEs, trait = "t1")
  expect_silent(report(geStab, outfile = tempfile(fileext = ".pdf")))
})

## varCov

test_that("varCov report functions correctly", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  skip_on_ci()
  geVarCov <- gxeVarCov(TD = BLUEs, trait = "t1")
  expect_silent(report(geVarCov, outfile = tempfile(fileext = ".pdf")))
})

