context("Summaries")

testTD <- createTD(data = testData, genotype = "seed", trial = "field",
                   rowId = "Y", colId = "X", rowCoord = "Y", colCoord = "X")
modelSp <- STRunModel(testTD, design = "rowcol", traits = c("t1", "t2"))
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

test_that("QTLDet summary produces correct output", {
  QTLDetF2_1 <- QTLDetect(testF2, trait = "phenotype")
  QTLDetF2_2 <- QTLDetect(testF2, trait = "phenotype", thrType = "fixed",
                          thrFixed = 5)
  sumQTLDetF2_1 <- capture.output(summary(QTLDetF2_1))
  sumQTLDetF2_2 <- capture.output(summary(QTLDetF2_2))
  expect_true("Peaks" %in% sumQTLDetF2_1)
  expect_true("No peaks detected" %in% sumQTLDetF2_2)
})
