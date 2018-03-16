context("utils")

test_that("function tryCatchExt catches errors", {
  catchErr <- tryCatchExt(stop("testErr"))
  expect_null(catchErr$value)
  expect_null(catchErr$warning)
  expect_is(catchErr$error, "error")
  expect_identical(catchErr$error$message, "testErr")
})

test_that("function tryCatchExt catches warnings", {
  catchWarn <- tryCatchExt(warning("testWng"))
  expect_identical(catchWarn$value, "testWng")
  expect_is(catchWarn$warning, "warning")
  expect_identical(catchWarn$warning$message, "testWng")
  expect_null(catchWarn$error)
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
  expect_is(catchWarnVal$warning, "warning")
  expect_identical(catchWarnVal$warning$message, "testWng")
  expect_null(catchWarnVal$error)
  catchWarnErr <- tryCatchExt({warning("testWng"); stop("testErr")})
  expect_null(catchWarnErr$value)
  expect_is(catchWarnErr$warning, "warning")
  expect_identical(catchWarnErr$warning$message, "testWng")
  expect_is(catchWarnErr$error, "error")
  expect_identical(catchWarnErr$error$message, "testErr")
})
