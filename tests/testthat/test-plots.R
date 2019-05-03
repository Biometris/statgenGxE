context("Plots")

## Testing the exact plot output is difficult but since also the ggplot
## objects on which the plots are based are invisibly returned at least some
## checking can be done.

geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
p0_1 <- plot(geAmmi, output = FALSE)
p0_2 <- plot(geAmmi, plotType = "AMMI2", output = FALSE)
test_that("AMMI plot gives correct output types", {
  expect_error(plot(geAmmi, plotType = "AMMI2", primAxis = "3"),
               "string starting with PC")
  expect_error(plot(geAmmi, plotType = "AMMI2", primAxis = "PC2"),
               "primAxis should differ from secAxis")
  expect_error(plot(geAmmi, plotType = "AMMI2", primAxis = "PC3"),
               "run with 2 principal components")
  expect_error(plot(geAmmi, plotType = "AMMI2", secAxis = "3"),
               "string starting with PC")
  expect_error(plot(geAmmi, plotType = "AMMI2", secAxis = "PC1"),
               "primAxis should differ from secAxis")
  expect_error(plot(geAmmi, plotType = "AMMI2", secAxis = "PC3"),
               "run with 2 principal components")
  expect_is(p0_1, "ggplot")
  expect_is(p0_2, "ggplot")
})

test_that("AMMI plot plotType options function properly", {
  geGGE <- gxeAmmi(TD = TDMaize, trait = "yld", GGE = TRUE)
  p1a <- plot(geAmmi, output = FALSE)
  p1b <- plot(geAmmi, plotType = "GGE1", output = FALSE)
  p2 <- plot(geGGE, output = FALSE)
  expect_equal(p1a, p1b)
  expect_equal(p2$labels$title, "GGE biplot for yld ")
})

test_that("AMMI plot scale option functions properly", {
  ## Only relevant for AMMI2 plots.
  p1_2 <- plot(geAmmi, plotType = "AMMI2", scale = 0, output = FALSE)
  p2_2 <- plot(geAmmi, plotType = "AMMI2", scale = 1, output = FALSE)
  p3_2 <- plot(geAmmi, plotType = "AMMI2", scale = 0.75, output = FALSE)
  expect_equal(p1_2$labels$title, "AMMI2 biplot for yld (genotype scaling) ")
  expect_equal(p2_2$labels$title, "AMMI2 biplot for yld (environment scaling) ")
  expect_equal(p3_2$labels$title, "AMMI2 biplot for yld (51%) ")
})

test_that("AMMI plot plotGeno functions properly", {
  p1_1 <- plot(geAmmi, plotGeno = FALSE, output = FALSE)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", plotGeno = FALSE, output = FALSE)
  ## Difference with default plot p0 should be the missing GeomPoint layer.
  geoms0_1 <- sapply(p0_1$layers, function(x) class(x$geom)[1])
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms0_2 <- sapply(p0_2$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_equal(setdiff(geoms0_1, geoms1_1), "GeomPoint")
  expect_equal(setdiff(geoms0_2, geoms1_2), "GeomPoint")
})

test_that("AMMI plot sizeGeno functions properly", {
  p1_1 <- plot(geAmmi, sizeGeno = 5, output = FALSE)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", sizeGeno = 5, output = FALSE)
  ## Difference with default plot p0 should be the replaced GeomPoint layer
  ## by GeomText layer.
  geoms0_1 <- sapply(p0_1$layers, function(x) class(x$geom)[1])
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms0_2 <- sapply(p0_2$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_equal(geoms1_1[geoms0_1 == "GeomPoint"], "GeomText")
  expect_equal(p1_1$layers[geoms0_1 == "GeomPoint"][[1]]$aes_params$size, 5)
  expect_equal(geoms1_2[geoms0_2 == "GeomPoint"], "GeomText")
  expect_equal(p1_2$layers[geoms0_2 == "GeomPoint"][[1]]$aes_params$size, 5)
})

test_that("AMMI plot plotEnv functions properly", {
  p1_1 <- plot(geAmmi, plotEnv = FALSE, output = FALSE)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", plotEnv = FALSE, output = FALSE)
  ## Difference with default plot p0 should be the missing GeomText layer.
  ## For AMMI2 also the arrows, GeomSegment, should be missing.
  geoms0_1 <- sapply(p0_1$layers, function(x) class(x$geom)[1])
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms0_2 <- sapply(p0_2$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_equal(setdiff(geoms0_1, geoms1_1), "GeomText")
  expect_setequal(setdiff(geoms0_2, geoms1_2), c("GeomText", "GeomSegment"))
})

test_that("AMMI plot sizeEnv functions properly", {
  p1_1 <- plot(geAmmi, sizeEnv = 5, output = FALSE)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", sizeEnv = 5, output = FALSE)
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_equal(p1_1$layers[geoms1_1 == "GeomText"][[1]]$aes_params$size, 5)
  expect_equal(p1_2$layers[geoms1_2 == "GeomText"][[1]]$aes_params$size, 5)
})

test_that("AMMI plot envFactor functions properly", {
  p1_1 <- plot(geAmmi, envFactor = 5, output = FALSE)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", envFactor = 5, output = FALSE)
  ## Coordinate limits should blow up by a factor 5.
  ## Lower x-limit is strongly dependent on genoscores so not blown up as much.
  expect_equal(5 * p0_1$plot_env$p$coordinates$limits$x[2],
               p1_1$plot_env$p$coordinates$limits$x[2])
  expect_equal(5 * p0_1$plot_env$p$coordinates$limits$y,
               p1_1$plot_env$p$coordinates$limits$y)
  ## Coordinate limits should blow up by a factor 5.
  expect_equal(5 * p0_2$plot_env$p$coordinates$limits$x,
               p1_2$plot_env$p$coordinates$limits$x)
  expect_equal(5 * p0_2$plot_env$p$coordinates$limits$y,
               p1_2$plot_env$p$coordinates$limits$y)
})

test_that("AMMI plot colEnv functions properly", {
  p1_1 <- plot(geAmmi, colEnv = "green", output = FALSE)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", colEnv = "green", output = FALSE)
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_equal(p1_1$layers[geoms1_1 == "GeomText"][[1]]$aes_params$colour, "green")
  expect_equal(p1_2$layers[geoms1_2 == "GeomText"][[1]]$aes_params$colour, "green")
  expect_equal(p1_2$layers[geoms1_2 == "GeomSegment"][[1]]$aes_params$colour, "green")
})

test_that("AMMI plot plotConvHull functions properly", {
  ## plotConvHull should be ignored for AMMI1.
  expect_equal(p0_1, plot(geAmmi, plotConvHull = TRUE, output = FALSE))
  ## For AMMI2 there should be two extra layers.
  p1_2 <- plot(geAmmi, plotType = "AMMI2", plotConvHull = TRUE, output = FALSE)
  geoms0_2 <- sapply(p0_2$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_setequal(geoms1_2[-match(geoms0_2, geoms1_2)],
                  c("GeomPolygon", "GeomSegment"))
})

testTDYear <- createTD(data = testDataYear, genotype = "seed",
                       trial = "field", repId = "rep",
                       subBlock = "block", rowId = "Y", colId = "X",
                       rowCoord = "Y", colCoord = "X")
modelSp <- fitTD(testTDYear, design = "rowcol", traits = c("t1", "t2"))
BLUEsYear <- SSAtoTD(modelSp, what = "BLUEs", keep = "year")
geAmmiYear <- gxeAmmi(BLUEsYear, trait = "t1", byYear = TRUE)
test_that("AMMI plot gives correct output types when byYear = TRUE", {
  p1 <- plot(geAmmiYear, output = FALSE)
  expect_is(p1, "list")
  expect_length(p1, 2)
  expect_named(p1, c("1", "2"))
  expect_is(p1[[1]], "ggplot")
  expect_is(p1[[2]], "ggplot")
  p2 <- plot(geAmmiYear, plotType = "AMMI2", output = FALSE)
  expect_is(p2, "list")
  expect_length(p2, 2)
  expect_named(p2, c("1", "2"))
  expect_is(p2[[1]], "ggplot")
  expect_is(p2[[2]], "ggplot")
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

SSA <- fitTD(TD = TDHeat05, design = "res.rowcol", traits = "yield")
test_that("SSA base plot gives correct output types", {
  p1 <- plot(SSA, traits = "yield", output = FALSE)
  expect_is(p1, "list")
  expect_length(p1, 1)
  expect_is(p1[[1]], "list")
  expect_length(p1[[1]], 1)
  expect_is(p1[[1]][[1]], "list")
  expect_length(p1[[1]][[1]], 4)
  lapply(X = p1[[1]][[1]], FUN = expect_is, "ggplot")
})

test_that("SSA spatial plot gives correct output types", {
  p1 <- plot(SSA, plotType = "spatial", traits = "yield", output = FALSE)
  expect_is(p1, "list")
  expect_length(p1, 1)
  expect_is(p1[[1]], "list")
  expect_length(p1[[1]], 1)
  expect_is(p1[[1]][[1]], "list")
  expect_length(p1[[1]][[1]], 6)
  lapply(X = p1[[1]][[1]], FUN = expect_is, "ggplot")
})

test_that("option what in SSA plot functions properly", {
  p1 <- plot(SSA, what = "random", output = FALSE)
  p2 <- plot(SSA, plotType = "spatial", what = "random", output = FALSE)
  expect_is(p1, "list")
  expect_equal(p2[[1]][[1]][[5]]$labels$title, "Genotypic BLUPs")
  expect_equal(p2[[1]][[1]][[6]]$labels$x, "Genotypic BLUPs")
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

test_that("varComp plot gives correct output types", {
  geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
  p <- plot(geVarComp, output = FALSE)
  geoms <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_is(p, "ggplot")
  expect_setequal(geoms, c("GeomTile", "GeomText"))
})

## melting data in the plot function caused an error when trials have a
## numerical value. This should not be the case.
test_that("varComp plot gives correct output types when trials are numerical", {
  TDMaize2 <- TDMaize
  for (trial in seq_along(TDMaize2)) {
    levels(TDMaize2[[trial]][["trial"]]) <- 1:8
  }
  geVarComp <- gxeVarComp(TD = TDMaize2, trait = "yld")
  expect_silent(p <- plot(geVarComp, output = FALSE))
})
