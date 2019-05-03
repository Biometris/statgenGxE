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
