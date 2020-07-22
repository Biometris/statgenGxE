context("Plots")

## Testing the exact plot output is difficult but since also the ggplot
## objects on which the plots are based are invisibly returned at least some
## checking can be done.

## AMMI

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
  p1a <- plot(geAmmi, plotType = "AMMI2")
  p1b <- plot(geAmmi, plotType = "GGE2")
  p2 <- plot(geGGE)
  expect_equal(p1a, p1b)
  expect_equal(p2$labels$title, "GGE biplot for t1 (environment scaling) ")
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
  geAmmi1 <- geAmmi
  ## Add missing values.
  geAmmi1$dat[geAmmi1$dat[["regime"]] == "W", "regime"] <- NA
  expect_error(plot(geAmmi, colorEnvBy = 1),
               "NULL or a character vector")
  expect_error(plot(geAmmi, colorEnvBy = "col"),
               "col has to be a column")
  expect_error(plot(geAmmi, colorEnvBy = "family"),
               "exactly one value per environment")
  expect_error(plot(geAmmi, colorEnvBy = "regime", colEnv = "blue"))
  expect_error(plot(geAmmi1, colorEnvBy = "regime"),
               "Missing values in regime")

  ## AMMI1
  p1_1 <- plot(geAmmi, colorEnvBy = "regime")
  p1_2 <- plot(geAmmi, colorEnvBy = "regime", colEnv = c("green", "blue"))
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(as.character(dat1_1[[".color"]]),
               c("#AA0DFE", "#AA0DFE", "#3283FE"))
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
               c("#AA0DFE", "#AA0DFE", "#3283FE"))
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomText"][[1]]$data
  expect_equal(as.character(dat1_2[[".color"]]),
               c("green", "green", "blue"))
})

test_that("AMMI plot colorGenoBy functions properly", {
  geAmmi1 <- geAmmi
  ## Add missing values.
  geAmmi1$dat[geAmmi1$dat[["family"]] == "F1", "family"] <- NA
  expect_error(plot(geAmmi, colorGenoBy = 1),
               "NULL or a character vector")
  expect_error(plot(geAmmi, colorGenoBy = "col"),
               "col has to be a column")
  expect_error(plot(geAmmi, colorGenoBy = "regime"),
               "exactly one value per genotype")
  expect_error(plot(geAmmi, colorGenoBy = "family", colGeno = "blue"))
  expect_error(plot(geAmmi1, colorGenoBy = "family"),
               "Missing values in family")

  ## AMMI1
  p1_1 <- plot(geAmmi, colorGenoBy = "family")
  p1_2 <- plot(geAmmi, colorGenoBy = "family",
               colGeno = c("green", "blue", "red"))
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomPoint"][[1]]$data
  expect_equal(as.character(dat1_1[[".color"]]),
               c("#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77",
                 "#D95F02", "#D95F02", "#D95F02", "#D95F02", "#D95F02",
                 "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#7570B3"))
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomPoint"][[1]]$data
  expect_equal(as.character(dat1_2[[".color"]]),
               c("green", "green", "green", "green", "green", "blue", "blue",
                 "blue", "blue", "blue", "red", "red", "red", "red", "red"))

  ## AMMI2
  p1_1 <- plot(geAmmi, plotType = "AMMI2", colorGenoBy = "family")
  p1_2 <- plot(geAmmi, plotType = "AMMI2", colorGenoBy = "family",
               colGeno = c("green", "blue", "red"))
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  dat1_1 <- p1_1$layers[geoms1_1 == "GeomPoint"][[1]]$data
  expect_equal(as.character(dat1_1[[".color"]]),
               c("#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77",
                 "#D95F02", "#D95F02", "#D95F02", "#D95F02", "#D95F02",
                 "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#7570B3"))
  dat1_2 <- p1_2$layers[geoms1_2 == "GeomPoint"][[1]]$data
  expect_equal(as.character(dat1_2[[".color"]]),
               c("green", "green", "green", "green", "green", "blue", "blue",
                 "blue", "blue", "blue", "red", "red", "red", "red", "red"))
})

test_that("colorEnvBy combined with colorGenoBy functions properly", {
  ## AMMI1
  p1_1 <- plot(geAmmi, colorGenoBy = "family", colorEnvBy = "regime")
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  datG1_1 <- p1_1$layers[geoms1_1 == "GeomPoint"][[1]]$data
  datE1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(as.character(datG1_1[[".color"]]),
               c("#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77",
                 "#D95F02", "#D95F02", "#D95F02", "#D95F02", "#D95F02",
                 "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#7570B3"))
  expect_equal(as.character(datE1_1[[".color"]]),
               c("#AA0DFE", "#AA0DFE", "#3283FE"))

  ## AMMI2
  p1_1 <- plot(geAmmi, plotType = "AMMI2", colorGenoBy = "family",
               colorEnvBy = "regime")
  geoms1_1 <- sapply(p1_1$layers, function(x) class(x$geom)[1])
  datG1_1 <- p1_1$layers[geoms1_1 == "GeomPoint"][[1]]$data
  datE1_1 <- p1_1$layers[geoms1_1 == "GeomText"][[1]]$data
  expect_equal(as.character(datG1_1[[".color"]]),
               c("#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77",
                 "#D95F02", "#D95F02", "#D95F02", "#D95F02", "#D95F02",
                 "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#7570B3"))
  expect_equal(as.character(datE1_1[[".color"]]),
               c("#AA0DFE", "#AA0DFE", "#3283FE"))
})

test_that("AMMI plot plotConvHull functions properly", {
  ## plotConvHull should be ignored for AMMI1.
  expect_equal(p0_1, plot(geAmmi, plotConvHull = TRUE))
  ## For AMMI2 there should be an extra layers.
  p1_2 <- plot(geAmmi, plotType = "AMMI2", plotConvHull = TRUE)
  geoms0_2 <- sapply(p0_2$layers, function(x) class(x$geom)[1])
  geoms1_2 <- sapply(p1_2$layers, function(x) class(x$geom)[1])
  expect_setequal(geoms1_2[-match(geoms0_2, geoms1_2)], "GeomPolygon")
})

test_that("GGE plot plotConvHull functions properly", {
  geGGE <- gxeGGE(TD = BLUEs, trait = "t1")
  ## For GGE2 there should be two extra layers.
  p_1 <- plot(geGGE, plotType = "GGE2")
  p_2 <- plot(geGGE, plotType = "GGE2", plotConvHull = TRUE)
  geoms_1 <- sapply(p_1$layers, function(x) class(x$geom)[1])
  geoms_2 <- sapply(p_2$layers, function(x) class(x$geom)[1])
  expect_setequal(geoms_2[-match(geoms_1, geoms_2)],
                  c("GeomPolygon", "GeomSegment"))
})

geAmmiYear <- gxeAmmi(BLUEsYear, trait = "t1", byYear = TRUE)
test_that("AMMI plot gives correct output types when byYear = TRUE", {
  geAmmiYear1 <- geAmmiYear2 <- geAmmiYear
  ## Add missing values.
  geAmmiYear1$dat[[1]][geAmmiYear1$dat[[1]][["regime"]] == "W", "regime"] <- NA
  geAmmiYear1$dat[[2]][geAmmiYear1$dat[[2]][["regime"]] == "W", "regime"] <- NA
  geAmmiYear2$dat[[1]][geAmmiYear2$dat[[1]][["family"]] == "F1", "family"] <- NA
  geAmmiYear2$dat[[2]][geAmmiYear2$dat[[2]][["family"]] == "F1", "family"] <- NA
  ## Year specific errors.
  expect_error(plot(geAmmiYear, colorGenoBy = "regime"),
               "exactly one value per genotype")
  expect_error(plot(geAmmiYear, colorEnvBy = "family"),
               "exactly one value per environment")
  expect_error(plot(geAmmiYear1, colorEnvBy = "regime"),
               "Missing values in regime")
  expect_error(plot(geAmmiYear2, colorGenoBy = "family"),
               "Missing values in family")
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

## Finlay Wilkinson

geFw <- gxeFw(TD = testTD, trait = "t1", maxIter = 30)

test_that("general check in FW plot function properly", {
  expect_error(plot(geFw, colorGenoBy = "col"), "col has to be a column")
})

test_that("FW plot gives correct output types", {
  p1 <- plot(geFw)
  p2 <- plot(geFw, plotType = "line")
  p3 <- plot(geFw, plotType = "trellis")
  p4 <- plot(geFw, plotType = "scatterFit")
  expect_is(p1, "list")
  expect_length(p1, 3)
  lapply(X = p1, FUN = expect_is, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
  expect_is(p4, "ggplot")
})

test_that("Option colorGenoBy in scatter plot functions correctly", {
  p1 <- plot(geFw, colorGenoBy = "family")
  expect_equal(p1[[1]]$labels$colour, "family")
  expect_equal(p1[[2]]$labels$colour, "family")
  expect_equal(p1[[3]]$labels$colour, "family")
})

test_that("Option colorGenoBy in line plot functions correctly", {
  p1 <- plot(geFw, plotType = "line", colorGenoBy = "family")
  expect_equal(p1$labels$colour, "family")
  ## With coloring plot should have a legend explicitly defined.
  expect_equal(p1$theme$legend.position, "right")
})

test_that("Option colorGenoBy in scatterFit plot functions correctly", {
  p1 <- plot(geFw, plotType = "scatterFit", colorGenoBy = "family")
  expect_equal(p1$labels$colour, "family")
})

test_that("option order in FW line plot functions properly", {
  p <- plot(geFw, plotType = "line", order = "descending")
  expect_equal(p$plot_env$xTrans, "reverse")
})

test_that("option genotypes in FW plot functions properly", {
  expect_error(plot(geFw, plotType = "trellis", genotypes = "g1"),
               "All genotypes should be in TD")
  p <- plot(geFw, plotType = "trellis", genotypes = paste0("G", 1:9))
  expect_equal(nlevels(p$data[["genotype"]]), 9)
})

## Stability

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

test_that("title functions correctly in stability plot", {
  geStab <- gxeStability(TD = testTD, trait = "t1")
  ## Actually just testing that it doesn't crash.
  ## Plots are returned as a list of plots,
  ## actual plotting, including title, is done by grid.arrange.
  expect_silent(plot(geStab, title = "Test"))
})

test_that("colorGenoBy functions correctly in stability plot", {
  geStab <- gxeStability(TD = testTD, trait = "t1")
  ## Actually just testing that it doesn't crash.
  ## Plots are returned as a list of plots,
  ## actual plotting, including legend, is done by grid.arrange
  expect_silent(plot(geStab, colorGenoBy = "family"))
})

## varCov

test_that("VarCov plot gives correct output types", {
  geVarCov <- gxeVarCov(TD = BLUEs, trait = "t1")
  p <- plot(geVarCov)
  geoms <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_is(p, "ggplot")
  expect_equal(geoms, "GeomTile")
})

## melting data in the plot function caused an error when trials have a
## numerical value. This should not be the case.
test_that("VarCov plot gives correct output types when trials are numerical", {
  for (trial in seq_along(BLUEs)) {
    levels(BLUEs[[trial]][["trial"]]) <- trial
  }
  geVarCov <- gxeVarCov(TD = BLUEs, trait = "t1")
  expect_silent(p <- plot(geVarCov))
})

## Mega environments.

geMegaEnv <- gxeMegaEnv(TD = BLUEs, trait = "t1")
test_that("megaEnv plot gives correct output types", {
  expect_warning(p <- plot(geMegaEnv),
                 "One should be cautious with the interpretation")
  expect_is(p, "list")
  expect_length(p, 1)
  expect_named(p, "pred")
  expect_is(p[[1]], "gtable")
  ## There should be 2 mega environments, so 2 x 2 panes in the plot layout.
  layout <- p[[1]]$layout
  expect_equal(nrow(layout[grepl(pattern  = "pane", x = layout[["name"]]), ]), 4)
})

test_that("option colorGenoBy in megaEnv plot functions correctly", {
  expect_warning(p0 <- plot(geMegaEnv))
  expect_warning(p <- plot(geMegaEnv, colorGenoBy = "family"),
                 "One should be cautious with the interpretation")
  ## New guide-box panel added.
  layout0 <- p0[[1]]$layout
  layout <- p[[1]]$layout
  expect_equal(setdiff(layout[["name"]], layout0[["name"]]), "guide-box")
})

## varComp

geVCLm <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4")
test_that("varComp plot gives correct output types", {
  p <- plot(geVCLm)
  expect_is(p, "ggplot")
})

test_that("option plotType in varComp plot functions correctly", {
  p0 <- plot(geVCLm)
  p1 <- plot(geVCLm, plotType = "percVar")
  expect_equal(p0$labels$x, "Square root of variance estimate")
  expect_equal(p1$labels$x, "Percentage of variance explained")
  expect_equal(p0$labels$title, "Standard deviations (general mean = 83)")
  expect_equal(p1$labels$title,
               "Percentage of variance explained (general mean = 83)")
})
