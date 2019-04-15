context("Report")

geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
test_that("Test that createReport function only accepts corretly named output", {
  expect_error(createReport(x = geAmmi, reportName = "ammiReport.Rnw",
                            outfile = tempfile(fileext = ".pd")),
               "invalid output filename")
  expect_error(createReport(x = geAmmi, reportName = "ammiReport.Rnw",
                            outfile = tempfile(fileext = "a b.pdf")),
               "outfile path cannot contain spaces")
  expect_silent(createReport(x = geAmmi, reportName = "ammiReport.Rnw",
                            outfile = tempfile(fileext = ".pdf")))
})
