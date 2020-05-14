context("Report")

geAmmi <- gxeAmmi(TD = BLUEs, trait = "t1")
test_that("Test that createReport function only accepts correctly named output", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  expect_error(report(geAmmi, outfile = tempfile(fileext = ".pd")),
               "Invalid output filename")
  expect_error(report(geAmmi, outfile = tempfile(fileext = "a b.pdf")),
               "outfile path cannot contain spaces")
  expect_silent(report(geAmmi, outfile = tempfile(fileext = ".pdf")))
})

