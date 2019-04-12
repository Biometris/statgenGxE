context("utils")

test_that("function tryCatchExt catches errors", {
  catchErr <- tryCatchExt(stop("testErr"))
  expect_null(catchErr$value)
  expect_null(catchErr$warning)
  expect_equal(catchErr$error, "testErr")
})

test_that("function tryCatchExt catches warnings", {
  catchWarn <- tryCatchExt(warning("testWng"))
  expect_equal(catchWarn$value, "testWng")
  expect_equal(catchWarn$warning, "testWng")
  expect_null(catchWarn$error)
  catch2Warn <- tryCatchExt({warning("testWng"); warning("testWng2")})
  expect_equal(catch2Warn$value, "testWng2")
  expect_equal(catch2Warn$warning, c("testWng", "testWng2"))
})

test_that("function tryCatchExt returns values", {
  catchVal <- tryCatchExt(1)
  expect_equal(catchVal$value, 1)
  expect_null(catchVal$warning)
  expect_null(catchVal$error)
})

test_that("function tryCatchExt returns combinations of outputs", {
  catchWarnVal <- tryCatchExt({warning("testWng"); 1})
  expect_equal(catchWarnVal$value, 1)
  expect_equal(catchWarnVal$warning, "testWng")
  expect_null(catchWarnVal$error)
  catchWarnErr <- tryCatchExt({warning("testWng"); stop("testErr")})
  expect_null(catchWarnErr$value)
  expect_equal(catchWarnErr$warning, "testWng")
  expect_equal(catchWarnErr$error, "testErr")
})

test_that("function supprWarn functions properly", {
  expect_warning(supprWarn(sqrt(-1), "testMsg"), "NaNs produced")
  expect_silent(supprWarn(sqrt(-1), "NaNs produced"))
})

test_that("function wrnToErr functions properly", {
  mod <- list(error = "testErr",
              warning = c("testWrn", "Abnormal termination"))
  modOut <- wrnToErr(mod)
  expect_equal(modOut$error, c("testErr", "Abnormal termination"))
  expect_equal(modOut$warning, "testWrn")
})

trDat <- data.frame(rowCoord = rep(1:4, each = 4),
                    colCoord = rep(1:4, times = 4),
                    repId = rep(1:2, each = 8))
test_that("calcPlotBorders functions properly", {
  bord <- calcPlotBorders(trDat = trDat, bordVar = "repId")
  expect_is(bord, "list")
  expect_named(bord, c("horW", "vertW"))
  expect_is(bord$horW, "data.frame")
  expect_is(bord$vertW, "data.frame")
  expect_equal(bord$horW$x, rep(1:4, each = 3))
  expect_equal(bord$horW$y, rep(c(1, 3, 5), times = 4))
  expect_equal(bord$vertW$x, rep(c(1, 5), times = 4))
  expect_equal(bord$vertW$y, rep(1:4, each = 2))
})

test_that("calcPlotBorders functions properly with NA", {
  trDat <- trDat[-1, ]
  bord <- calcPlotBorders(trDat = trDat, bordVar = "repId")
  expect_equal(bord$horW$x, c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4))
  expect_equal(bord$horW$y, c(3, 5, 1, 3, 5, 1, 3, 5, 1, 3, 5))
  expect_equal(bord$vertW$x, c(5, 1, 5, 1, 5, 1, 5))
  expect_equal(bord$vertW$y, c(1, 2, 2, 3, 3, 4, 4))
})

test_that("mapData functions properly", {
  mapDat <- mapData(xLim = c(0, 5), yLim = c(50, 53))
  expect_is(mapDat, "data.frame")
  expect_true(all(hasName(mapDat, c("lat", "long"))))
  expect_equal(sum(is.na(mapDat$lat)), 0)
  expect_equal(sum(is.na(mapDat$long)), 0)
  expect_equal(unique(mapDat$region),
                      c("Belgium", "France", "UK", "Netherlands"))
})
