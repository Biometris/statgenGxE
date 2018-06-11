context("multiQTLFit")

QTLDetF2 <- QTLDetect(testF2, trait = "phenotype", thrType = "fixed",
                      thrFixed = 1.5, window = 2)
mqf <- multiQTLFit(QTLDetF2)
test_that("multiQTLFit returns proper output format", {
  expect_is(mqf, "multiQTL")
  expect_is(mqf$qtl, "fitqtl")
  expect_equal(mqf$QTLDet, QTLDetF2)
})

test_that("multiQTLFit returns proper output values", {
  expect_equivalent(unlist(mqf$qtl$result.full[1, ]),
                    c(2, 261.981621283948, 130.990810641974, 2.86346437478481,
                      23.1822778541954, 0.00136941671715973,
                      0.00203396505756059))
  expect_equal(names(mqf$qtl$ests$ests), c("Intercept", "D1M318a", "D1M318d"))
  expect_equivalent(mqf$qtl$ests$ests,
                    c(24.1128838934645, 3.81180377175034, -0.34710296634032))
  expect_equal(colnames(mqf$qtl$ests$covar), c("Intercept", "1@9.8a", "1@9.8d"))
  expect_equivalent(mqf$qtl$ests$covar,
                    c(0.455820146628824, 0.314854311126121, -0.238279871090813,
                      0.314854311126121, 1.18364334326509, -0.65054406275895,
                      -0.238279871090813, -0.65054406275895, 1.95857006243086))
  expect_equal(mqf$qtl$lod, 2.86346437478481)
})

test_that("option selection functions properly", {
  mqfSel <- multiQTLFit(QTLDetF2, selection = "none")
  expect_equal(rownames(mqfSel$qtl$result.drop),
               c("D1M212", "D1M318", "D2M241", "D2M336", "DXM66"))
  expect_equal(names(mqfSel$qtl$ests$ests),
               c("Intercept", "D1M212.2", "D1M212.3", "D1M318.2", "D1M318.3",
                 "D2M241.2", "D2M241.3", "D2M336.2", "D2M336.3", "DXM66.2",
                 "DXM66.3", "DXM66.4"))
  expect_equivalent(mqfSel$qtl$ests$ests,
                    c(17.6873150122593, 5.59009305825471, -1.91468170015809,
                      -2.80487463939015, 8.33269214212498, 2.23120858065274,
                      1.71736319942582, -0.707650205305372, 1.64968633500598,
                      3.04748476414923, 0.263195708577232, 3.21169324590292))
})
