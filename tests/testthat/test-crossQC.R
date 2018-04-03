context("crossQC")

test_that("quality control on cross objects produces new cross object", {
  expect_error(QTLMapQC(testF2), "You're attempting to drop")
  ## Warning caused by low number of individuals in check segDistortion.
  expect_warning(QTLMapQC(testF2, missMrk = 0.2),
                 "Chi-squared approximation may be incorrect")
  expect_is(QTLMapQC(testF2, missMrk = 0.2, missInd = 0.2), "cross")
})

test_that("option missMrk functions properly", {
  expect_equivalent(qtl::nmar(QTLMapQC(testF2, missMrk = 0, missInd = 0)),
                    c(7, 5, 3))
  expect_equivalent(qtl::nmar(QTLMapQC(testF2, missMrk = 0.2, missInd = 0)),
                    c(4, 5, 2))
})

test_that("option missInd functions properly", {
  expect_equal(qtl::nind(QTLMapQC(testF2, missMrk = 0, missInd = 0)), 41)
  expect_equal(qtl::nind(QTLMapQC(testF2, missMrk = 0, missInd = 0.2)), 30)
})

test_that("option removeDuplicates functions properly", {
  expect_equivalent(qtl::nmar(QTLMapQC(testF2, missMrk = 0, missInd = 0,
                                       removeDuplicates = FALSE)), c(8, 5, 3))
})

test_that("option segDistortion functions properly", {
  expect_equivalent(qtl::nmar(QTLMapQC(testF2, missMrk = 0.2,
                                       segDistortion = 0)), c(4, 5, 2))
  expect_equivalent(qtl::nmar(QTLMapQC(testF2, missMrk = 0.2, missInd = 0,
                                       segDistortion = 0.2)), c(3, 5, 2))
})

test_that("option recombination functions properly", {
  testF2Sw <- qtl::switchAlleles(testF2, "D1M430")
  expect_equivalent(qtl::nmar(QTLMapQC(testF2Sw, missMrk = 0, missInd = 0,
                                       recombination = 0)), c(7, 5, 3))
  expect_equivalent(qtl::nmar(QTLMapQC(testF2Sw, missMrk = 0, missInd = 0,
                                       recombination = 3)), c(6, 5, 3))
})

test_that("options reestimateMap functions properly", {
  testF2Est <- QTLMapQC(testF2, missMrk = 0, missInd = 0, reestimateMap = TRUE)
  expect_equivalent(qtl::nmar(testF2Est), c(7, 5, 3))
  expect_equal(qtl::nind(testF2Est), 41)
  expect_equal(qtl::map2table(qtl::pull.map(testF2Est))$pos,
               c(0, 11.2884499432804, 13.1079823171321, 56.240356992947,
                 73.5189709312085, 92.4970905547784, 104.828719044552,
                 0, 29.8623166499185, 48.2035470096869, 57.802867410675,
                 74.891205717183, 0, 9.81381920438241, 21.4203515860796))
})

test_that("option crossover functions properly", {
  expect_equal(qtl::nind(QTLMapQC(testF2, missMrk = 0, missInd = 0,
                                  crossover = 0)), 50)
  expect_equal(qtl::nind(QTLMapQC(testF2, missMrk = 0, missInd = 0,
                                  crossover = 0.3)), 46)
})

