context("Report")

geAmmi <- gxeAmmi(TD = BLUEs, trait = "t1")
test_that("Test that createReport function only accepts corretly named output", {
  ## Reporting requires pdflatex which isn't available on cran.
  skip_on_cran()
  expect_error(createReport(x = geAmmi, reportName = "ammiReport.Rnw",
                            outfile = tempfile(fileext = ".pd")),
               "invalid output filename")
  expect_error(createReport(x = geAmmi, reportName = "ammiReport.Rnw",
                            outfile = tempfile(fileext = "a b.pdf")),
               "outfile path cannot contain spaces")
  expect_silent(createReport(x = geAmmi, reportName = "ammiReport.Rnw",
                            outfile = tempfile(fileext = ".pdf")))
})
