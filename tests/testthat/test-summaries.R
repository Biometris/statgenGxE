context("Summaries")

test_that("AMMI summary produces correct output", {
  geAmmi <- gxeAmmi(TD = BLUEs, trait = "t1")
  geGGE <- gxeGGE(TD = BLUEs, trait = "t1")
  sumAmmi <- capture.output(summary(geAmmi))
  sumGGE <- capture.output(summary(geGGE))
  sumAmmi2 <- capture.output(summary(geAmmi, printGenoScores = TRUE))
  expect_true(all(c("Principal components ", "Anova ",
                    "Environment scores ") %in% sumAmmi))
  expect_false("Genotypic scores " %in% sumAmmi)
  expect_true(all(c("Principal components ", "Environment scores ") %in% sumAmmi))
  expect_false(all(c("Anova ", "Genotypic scores ") %in% sumGGE))
  expect_true(all(c("Principal components ", "Anova ",
                    "Environment scores ", "Genotypic scores ") %in% sumAmmi2))
})

test_that("AMMI summary produces correct output per year", {
  geAmmiYear <- gxeAmmi(BLUEsYear, trait = "t1", byYear = TRUE)
  sumAmmiYear <- capture.output(summary(geAmmiYear, printGenoScores = TRUE))
  ## Checking that output is printed for both years.
  expect_length(grep("Standard deviation", sumAmmiYear), 2)
  expect_length(grep("Interactions", sumAmmiYear), 2)
  expect_length(grep("E1", sumAmmiYear), 1)
  expect_length(grep("E4", sumAmmiYear), 1)
  expect_length(grep("G2", sumAmmiYear), 2)
})

test_that("FW summary produces correct output", {
  ## t1 doesn't converge. Use t2 instead.
  geFW <- gxeFw(TD = BLUEs, trait = "t2")
  sumFW <- capture.output(summary(geFW))
  geFW2 <- gxeFw(TD = BLUEs, trait = "t2", sorted = "ascending")
  sumFW2 <- capture.output(summary(geFW2))
  geFW3 <- gxeFw(TD = BLUEs, trait = "t2", sorted = "none")
  sumFW3 <- capture.output(summary(geFW3))
  expect_true(all(c("Environmental effects ", "Anova ",
                    "Most sensitive genotypes") %in% sumFW))
  expect_true("Least sensitive genotypes" %in% sumFW2)
  expect_true("First five genotypes" %in% sumFW3)
})

test_that("varCov summary produces correct output", {
  geVC <- gxeVarCov(TD = BLUEs, trait = "t1")
  sumVC <- capture.output(summary(geVC))
  expect_true("Best model: cs, based on BIC." %in% sumVC)
})

test_that("Stability summary produces correct output", {
  geStab <- gxeStability(TD = BLUEs, trait = "t1")
  sumStab <- capture.output(summary(geStab))
  sumStab2 <- capture.output(summary(geStab, pctGeno = 20))
  expect_true(all(c("Cultivar-superiority measure (Top 10 % genotypes)",
                    "Static stability (Top 10 % genotypes)",
                    "Wricke's ecovalence (Top 10 % genotypes)") %in% sumStab))
  expect_equal(length(sumStab2), length(sumStab) + 3)
})

test_that("megaEnv summary produces correct output", {
  geMegaEnv <- gxeMegaEnv(TD = BLUEs, trait = "t1")
  sumMegaEnv <- capture.output(summary(geMegaEnv))
  expect_true("Mega environments based on t1" %in% sumMegaEnv)
  expect_true(" Mega_factor Trial Winning_genotype AMMI_estimates" %in% sumMegaEnv)
})

test_that("varComp summary produces correct output", {
  geVCLm <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4")
  sumVCLm <- capture.output(summary(geVCLm))
  expect_true("t1 ~ trial + (1 | genotype) " %in% sumVCLm)
  expect_true("Sources of variation" %in% sumVCLm)

  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
  sumVCAs <- capture.output(summary(geVCAs))
  expect_true("t1 ~ trial + (1 | genotype) " %in% sumVCAs)
  expect_true("Sources of variation" %in% sumVCAs)
})

