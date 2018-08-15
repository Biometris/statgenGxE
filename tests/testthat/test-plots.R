context("Plots")

## Testing the exact plot output is difficult but since also the ggplot
## objects on which the plots are based are invisibly returned at least some
## checking can be done.

test_that("AMMI plot gives correct output types", {
  geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
  p1 <- plot(geAmmi, output = FALSE)
  expect_error(plot(geAmmi, plotType = "AMMI2", secAxis = "3"),
               "string starting with PC")
  expect_error(plot(geAmmi, plotType = "AMMI2", secAxis = "PC1"),
               "Invalid value provided for secAxis")
  expect_error(plot(geAmmi, plotType = "AMMI2", secAxis = "PC3"),
               "run with 2 principal components")
  p2 <- plot(geAmmi, plotType = "AMMI2", output = FALSE)
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
})

test_that("FW plot gives correct output types", {
  geFw <- geFW <- gxeFw(TD = TDMaize, trait = "yld")
  p1 <- plot(geFw, output = FALSE)
  p2 <- plot(geFw, plotType = "line", output = FALSE)
  p3 <- plot(geFw, plotType = "trellis", output = FALSE)
  expect_is(p1, "list")
  expect_length(p1, 3)
  lapply(X = p1, FUN = expect_is, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
})

test_that("SSA plot gives correct output types", {
  SSA <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield")
  p1 <- plot(SSA, output = FALSE)
  p2 <- plot(SSA, plotType = "spatial", output = FALSE)
  expect_is(p1, "list")
  expect_length(p1, 4)
  lapply(X = p1, FUN = expect_is, "ggplot")
  expect_is(p2, "list")
  expect_length(p2, 6)
  lapply(X = p2, FUN = expect_is, "ggplot")
})

test_that("stability plot gives correct output types", {
  geStab <- gxeStability(TD = TDMaize, trait = "yld")
  p1 <- plot(geStab, output = FALSE)
  expect_is(p1, "list")
  expect_length(p1, 4)
  lapply(X = p1, FUN = expect_is, "ggplot")
  geStab2 <- gxeStability(TD = TDMaize, trait = "yld", method = "superiority")
  p2 <- plot(geStab2, output = FALSE)
  expect_length(p2, 1)
})

test_that("QTLDet plot gives correct output types", {
  QTLDetF2 <- QTLDetect(testF2, trait = "phenotype", thrType = "fixed",
                        thrFixed = 1.5, window = 2)
  p <- plot(QTLDetF2, output = FALSE)
  expect_is(p, "ggplot")
})

test_that("multiQTL plot gives correct output types", {
  QTLDetF2 <- QTLDetect(testF2, trait = "phenotype", thrType = "fixed",
                        thrFixed = 1.5, window = 2)
  mqf <- multiQTLFit(QTLDetF2)
  p <- plot(mqf, output = FALSE)
  expect_is(p, "ggplot")
})


