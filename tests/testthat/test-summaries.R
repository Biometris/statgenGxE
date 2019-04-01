context("Summaries")

test_that("AMMI summary produces correct output", {
  geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
  geGGE <- gxeAmmi(TD = TDMaize, trait = "yld", GGE = TRUE)
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
