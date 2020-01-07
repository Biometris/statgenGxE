context("Plots")

## Testing the exact plot output is difficult but since also the ggplot
## objects on which the plots are based are invisibly returned at least some
## checking can be done.

modelSp <- fitTD(testTD, design = "rowcol", traits = c("t1", "t2"))
BLUEs <- STAtoTD(modelSp, what = "BLUEs", keep = c("family", "regime"))
geAmmi <- gxeAmmi(BLUEs, trait = "t1")

test_that("general checks in ammi plot function properly", {
  expect_error(plot(geAmmi, scale = 2),
               "a single numerical value between 0 and 1")
  expect_error(plot(geAmmi, sizeGeno = -1),
               "a single numerical value greater than or equal to 0")
  expect_error(plot(geAmmi, sizeEnv = -1),
               "a single numerical value greater than or equal to 0")
  expect_error(plot(geAmmi, colGeno = 1),
               "NULL or a character vector")
  expect_error(plot(geAmmi, colEnv = 1),
               "NULL or a character vector")
  expect_error(plot(geAmmi, plotType = "AMMI2", primAxis = "PCC"),
               "Invalid value provided for primAxis")
  expect_error(plot(geAmmi, plotType = "AMMI2", secAxis = "PCC"),
               "Invalid value provided for secAxis")
})

p0_1 <- plot(geAmmi)
p0_2 <- plot(geAmmi, plotType = "AMMI2")
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
  geGGE <- gxeGGE(TD = BLUEs, trait = "t1")
  p1a <- plot(geAmmi)
  p1b <- plot(geAmmi, plotType = "GGE1")
  p2 <- plot(geGGE)
  expect_equal(p1a, p1b)
  expect_equal(p2$labels$title, "GGE biplot for t1 ")
})

test_that("AMMI plot scale option functions properly", {
  ## Only relevant for AMMI2 plots.
  p1_2 <- plot(geAmmi, plotType = "AMMI2", scale = 0)
  p2_2 <- plot(geAmmi, plotType = "AMMI2", scale = 1)
  p3_2 <- plot(geAmmi, plotType = "AMMI2", scale = 0.75)
  p4_2 <- plot(geAmmi, plotType = "AMMI2", scale = 0.5)
  expect_equal(p1_2$labels$title, "AMMI2 biplot for t1 (genotype scaling) ")
  expect_equal(p2_2$labels$title, "AMMI2 biplot for t1 (environment scaling) ")
  expect_equal(p3_2$labels$title, "AMMI2 biplot for t1 (100%) ")
  expect_equal(p4_2$labels$title, "AMMI2 biplot for t1 (symmetric scaling) ")
})

test_that("AMMI plot plotGeno functions properly", {
  p1_1 <- plot(geAmmi, plotGeno = FALSE)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", plotGeno = FALSE)
  ## Difference with default plot p0 should be the missing GeomPoint layer.
  geoms0_1 <- sapply(p0_1$layers, function(x) class(x$geom)[1])
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms0_2 <- sapply(p0_2$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_equal(setdiff(geoms0_1, geoms1_1), "GeomPoint")
  expect_equal(setdiff(geoms0_2, geoms1_2), "GeomPoint")
})

test_that("AMMI plot sizeGeno functions properly", {
  p1_1 <- plot(geAmmi, sizeGeno = 5)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", sizeGeno = 5)
  ## Difference with default plot p0 should be the replaced GeomPoint layer
  ## by GeomText layer.
  geoms0_1 <- sapply(p0_1$layers, function(x) class(x$geom)[1])
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms0_2 <- sapply(p0_2$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_equal(geoms1_1[geoms0_1 == "GeomPoint"], "GeomText")
  dat1_1 <- p1_1$layers[geoms0_1 == "GeomPoint"][[1]]$data
  expect_equal(unique(dat1_1[dat1_1[["type"]] == "geno", ".size"]), 5)
  expect_equal(geoms1_2[geoms0_2 == "GeomPoint"], "GeomText")
  dat1_2 <- p1_2$layers[geoms0_2 == "GeomPoint"][[1]]$data
  expect_equal(unique(dat1_1[dat1_1[["type"]] == "geno", ".size"]), 5)
})

test_that("AMMI plot plotEnv functions properly", {
  p1_1 <- plot(geAmmi, plotEnv = FALSE)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", plotEnv = FALSE)
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
  p1_1 <- plot(geAmmi, sizeEnv = 5)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", sizeEnv = 5)
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(unique(dat1_1[dat1_1[["type"]] == "env", ".size"]), 5)
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomText"][[1]]$data
  expect_equal(unique(dat1_2[dat1_2[["type"]] == "env", ".size"]), 5)
})

test_that("AMMI plot envFactor functions properly", {
  p1_1 <- plot(geAmmi, envFactor = 5)
  p1_2 <- plot(geAmmi, plotType = "AMMI2", envFactor = 5)
  ## Coordinate limits should blow up by a factor 5.
  ## x-limits is strongly dependent on genoscores so not blown up as much.
  expect_equal(5 * p0_1$plot_env$p$coordinates$limits$y,
               p1_1$plot_env$p$coordinates$limits$y)
  ## Coordinate limits should blow up by a factor 5.
  expect_equal(5 * p0_2$plot_env$p$coordinates$limits$x,
               p1_2$plot_env$p$coordinates$limits$x)
  expect_equal(5 * p0_2$plot_env$p$coordinates$limits$y,
               p1_2$plot_env$p$coordinates$limits$y)
})

test_that("AMMI plot colEnv functions properly", {
  p1_1 <- plot(geAmmi, colEnv = "green")
  p1_2 <- plot(geAmmi, plotType = "AMMI2", colEnv = "green")
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(as.character(unique(dat1_1[dat1_1[["type"]] == "env", ".color"])),
               "green")
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomText"][[1]]$data
  expect_equal(as.character(unique(dat1_2[dat1_2[["type"]] == "env", ".color"])),
               "green")
})

test_that("AMMI plot colorEnvBy functions properly", {
  expect_error(plot(geAmmi, colorEnvBy = 1),
               "NULL or a character vector")
  expect_error(plot(geAmmi, colorEnvBy = "col"),
               "col has to be a column")
  expect_error(plot(geAmmi, colorEnvBy = "family"),
               "exactly one value per environment")
  expect_error(plot(geAmmi, colorEnvBy = "regime", colEnv = "blue"))

  ## AMMI1
  p1_1 <- plot(geAmmi, colorEnvBy = "regime")
  p1_2 <- plot(geAmmi, colorEnvBy = "regime", colEnv = c("green", "blue"))
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(as.character(dat1_1[[".color"]]),
               c("#4C00FFFF", "#4C00FFFF", "#00E5FFFF"))
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomText"][[1]]$data
  expect_equal(as.character(dat1_2[[".color"]]),
               c("green", "green", "blue"))

  ## AMMI2
  p1_1 <- plot(geAmmi, plotType = "AMMI2", colorEnvBy = "regime")
  p1_2 <- plot(geAmmi, plotType = "AMMI2", colorEnvBy = "regime",
               colEnv = c("green", "blue"))
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(as.character(dat1_1[[".color"]]),
               c("#4C00FFFF", "#4C00FFFF", "#00E5FFFF"))
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomText"][[1]]$data
  expect_equal(as.character(dat1_2[[".color"]]),
               c("green", "green", "blue"))
})

test_that("AMMI plot colorGenoBy functions properly", {
  expect_error(plot(geAmmi, colorGenoBy = 1),
               "NULL or a character vector")
  expect_error(plot(geAmmi, colorGenoBy = "col"),
               "col has to be a column")
  expect_error(plot(geAmmi, colorGenoBy = "regime"),
               "exactly one value per genotype")
  expect_error(plot(geAmmi, colorGenoBy = "family", colGeno = "blue"))

  ## AMMI1
  p1_1 <- plot(geAmmi, colorGenoBy = "family")
  p1_2 <- plot(geAmmi, colorGenoBy = "family",
               colGeno = c("green", "blue", "red"))
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomPoint"][[1]]$data
  expect_equal(as.character(dat1_1[[".color"]]),
               c("#4C00FFFF", "#00FF4DFF", "#FFFF00FF", "#FFFF00FF",
                 "#FFFF00FF", "#FFFF00FF", "#FFFF00FF", "#4C00FFFF",
                 "#4C00FFFF", "#4C00FFFF", "#4C00FFFF", "#00FF4DFF",
                 "#00FF4DFF", "#00FF4DFF", "#00FF4DFF"))
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomPoint"][[1]]$data
  expect_equal(as.character(dat1_2[[".color"]]),
               c("green", "blue", "red", "red", "red", "red", "red", "green",
                 "green", "green", "green", "blue", "blue", "blue", "blue"))

  ## AMMI2
  p1_1 <- plot(geAmmi, plotType = "AMMI2", colorGenoBy = "family")
  p1_2 <- plot(geAmmi, plotType = "AMMI2", colorGenoBy = "family",
               colGeno = c("green", "blue", "red"))
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomPoint"][[1]]$data
  expect_equal(as.character(dat1_1[[".color"]]),
               c("#4C00FFFF", "#00FF4DFF", "#FFFF00FF", "#FFFF00FF",
                 "#FFFF00FF", "#FFFF00FF", "#FFFF00FF", "#4C00FFFF",
                 "#4C00FFFF", "#4C00FFFF", "#4C00FFFF", "#00FF4DFF",
                 "#00FF4DFF", "#00FF4DFF", "#00FF4DFF"))
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomPoint"][[1]]$data
  expect_equal(as.character(dat1_2[[".color"]]),
               c("green", "blue", "red", "red", "red", "red", "red", "green",
                 "green", "green", "green", "blue", "blue", "blue", "blue"))
})

test_that("colorEnvBy combined with colorGenoBy functions properly", {
  ## AMMI1
  p1_1 <- plot(geAmmi, colorGenoBy = "family", colorEnvBy = "regime")
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  datG1_1 <- p1_1$layers[geoms1_1 == "GeomPoint"][[1]]$data
  datE1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(as.character(datG1_1[[".color"]]),
               c("#4C00FFFF", "#00FF4DFF", "#FFFF00FF", "#FFFF00FF",
                 "#FFFF00FF", "#FFFF00FF", "#FFFF00FF", "#4C00FFFF",
                 "#4C00FFFF", "#4C00FFFF", "#4C00FFFF", "#00FF4DFF",
                 "#00FF4DFF", "#00FF4DFF", "#00FF4DFF"))
  expect_equal(as.character(datE1_1[[".color"]]),
               c("#4C00FFFF", "#4C00FFFF", "#00E5FFFF"))

  ## AMMI2
  p1_1 <- plot(geAmmi, plotType = "AMMI2", colorGenoBy = "family",
               colorEnvBy = "regime")
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  datG1_1 <- p1_1$layers[geoms1_1 == "GeomPoint"][[1]]$data
  datE1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(as.character(datG1_1[[".color"]]),
               c("#4C00FFFF", "#00FF4DFF", "#FFFF00FF", "#FFFF00FF",
                 "#FFFF00FF", "#FFFF00FF", "#FFFF00FF", "#4C00FFFF",
                 "#4C00FFFF", "#4C00FFFF", "#4C00FFFF", "#00FF4DFF",
                 "#00FF4DFF", "#00FF4DFF", "#00FF4DFF"))
  expect_equal(as.character(datE1_1[[".color"]]),
               c("#4C00FFFF", "#4C00FFFF", "#00E5FFFF"))
})

test_that("AMMI plot plotConvHull functions properly", {
  ## plotConvHull should be ignored for AMMI1.
  expect_equal(p0_1, plot(geAmmi, plotConvHull = TRUE))
  ## For AMMI2 there should be two extra layers.
  p1_2 <- plot(geAmmi, plotType = "AMMI2", plotConvHull = TRUE)
  geoms0_2 <- sapply(p0_2$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_setequal(geoms1_2[-match(geoms0_2, geoms1_2)],
                  c("GeomPolygon", "GeomSegment"))
})

modelSp <- fitTD(testTDYear, design = "rowcol", traits = c("t1", "t2"))
BLUEsYear <- STAtoTD(modelSp, what = "BLUEs",
                     keep = c("year", "family", "regime"))
geAmmiYear <- gxeAmmi(BLUEsYear, trait = "t1", byYear = TRUE)
test_that("AMMI plot gives correct output types when byYear = TRUE", {
  ## Year specific errors.
  expect_error(plot(geAmmiYear, colorGenoBy = "regime"),
               "exactly one value per genotype")
  expect_error(plot(geAmmiYear, colorEnvBy = "family"),
               "exactly one value per environment")
  expect_error(plot(geAmmiYear, plotType = "AMMI2", primAxis = "PC3"),
               "Highest number of principal components is 2")
  expect_error(plot(geAmmiYear, plotType = "AMMI2", secAxis = "PC3"),
               "Highest number of principal components is 2")

  p1 <- plot(geAmmiYear)
  expect_is(p1, "list")
  expect_length(p1, 2)
  expect_named(p1, c("1", "2"))
  expect_is(p1[[1]], "ggplot")
  expect_is(p1[[2]], "ggplot")
  p2 <- plot(geAmmiYear, plotType = "AMMI2")
  expect_is(p2, "list")
  expect_length(p2, 2)
  expect_named(p2, c("1", "2"))
  expect_is(p2[[1]], "ggplot")
  expect_is(p2[[2]], "ggplot")
})

geFw <- gxeFw(TD = testTD, trait = "t1", maxIter = 30)
test_that("FW plot gives correct output types", {
  p1 <- plot(geFw)
  p2 <- plot(geFw, plotType = "line")
  p3 <- plot(geFw, plotType = "trellis")
  expect_is(p1, "list")
  expect_length(p1, 3)
  lapply(X = p1, FUN = expect_is, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
})

test_that("option order in FW plot functions properly", {
  p <- plot(geFw, plotType = "line", order = "descending")
  expect_equal(p$plot_env$xTrans, "reverse" )
})

test_that("option genotypes in FW plot functions properly", {
  expect_error(plot(geFw, plotType = "trellis", genotypes = "g1"),
               "All genotypes should be in TD")
  p <- plot(geFw, plotType = "trellis", genotypes = paste0("G", 1:9))
  expect_equal(nlevels(p$data[["genotype"]]), 9)
})

test_that("stability plot gives correct output types", {
  geStab <- gxeStability(TD = testTD, trait = "t1")
  p1 <- plot(geStab)
  expect_is(p1, "list")
  expect_length(p1, 4)
  lapply(X = p1, FUN = expect_is, "ggplot")
  geStab2 <- gxeStability(TD = testTD, trait = "t1", method = "superiority")
  p2 <- plot(geStab2)
  expect_length(p2, 1)
})

test_that("title argument functions correctly in stability plot", {
  geStab <- gxeStability(TD = testTD, trait = "t1")
  ## Actually just testing that it doesn't crash.
  ## Plots are returned as a list of plots,
  ## actual plotting, including title, is done by grid.arrange.s
  expect_silent(plot(geStab, title = "Test"))
})

test_that("varComp plot gives correct output types", {
  geVarComp <- gxeVarComp(TD = testTD, trait = "t1")
  p <- plot(geVarComp)
  geoms <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_is(p, "ggplot")
  expect_setequal(geoms, c("GeomTile", "GeomText"))
})

## melting data in the plot function caused an error when trials have a
## numerical value. This should not be the case.
test_that("varComp plot gives correct output types when trials are numerical", {
  for (trial in seq_along(testTD)) {
    levels(testTD[[trial]][["trial"]]) <- 1:3
  }
  geVarComp <- gxeVarComp(TD = testTD, trait = "t1")
  expect_silent(p <- plot(geVarComp))
})
