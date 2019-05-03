context("Summaries")

testTD <- createTD(data = testData, genotype = "seed", trial = "field",
                   rowId = "Y", colId = "X", rowCoord = "Y", colCoord = "X")
modelSp <- fitTD(testTD, design = "rowcol", traits = c("t1", "t2"))
BLUEs <- SSAtoTD(modelSp, what = "BLUEs")

test_that("AMMI summary produces correct output", {
  geAmmi <- gxeAmmi(TD = BLUEs, trait = "t1")
  geGGE <- gxeAmmi(TD = BLUEs, trait = "t1", GGE = TRUE)
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

test_that("varComp summary produces correct output", {
  geVC <- gxeVarComp(TD = BLUEs, trait = "t1")
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