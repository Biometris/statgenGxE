context("gxeAmmi")

test_that("general checks in gxeAmmi function properly", {
  expect_error(gxeAmmi(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeAmmi(testTD, trait = c("t1", "t1")),
               "trait has to be a character string of length 1")
  expect_error(gxeAmmi(BLUEs, trait = "t5"),
               "t5 has to be a column in TD")
  expect_error(gxeAmmi(BLUEs, trials = "E4", trait = "t1"),
               "a character vector defining trials in TD")
  expect_error(gxeAmmi(BLUEs, trials = c("E1", "E2"), trait = "t1"),
               "contain at least 3 trials")
  expect_error(gxeAmmi(BLUEs, trait = "t1", byYear = TRUE),
               "year has to be a column in TD")
  expect_error(gxeAmmi(BLUEs, trait = "t1", useWt = TRUE),
               "wt has to be a column in TD")
  expect_error(gxeAmmi(BLUEs, trait = "t1", nPC = 0),
               "NULL or a single numerical value greater than or equal to 1")
})

test_that("check for proper GxE data function properly", {
  ## Duplicate first observation in the data for E1.
  BLUEs$E1 <- rbind(BLUEs$E1, BLUEs$E1[1, ])
  expect_error(gxeAmmi(BLUEs, trait = "t1"),
               "at most 1 value per trial per genotype")
  BLUEsYear$E1 <- rbind(BLUEsYear$E1, BLUEsYear$E1[1, ])
  expect_warning(gxeAmmi(BLUEsYear, trait = "t1", byYear = TRUE),
                 "More than 1 value per trial per genotype for 1")
})

geAmmi <- gxeAmmi(BLUEs, trait = "t1")
test_that("output is of the right class", {
  expect_is(geAmmi, "AMMI")
  expect_is(geAmmi$envScores, "matrix")
  expect_is(geAmmi$genoScores, "matrix")
  expect_is(geAmmi$importance, "data.frame")
  expect_is(geAmmi$anova, "data.frame")
  expect_is(geAmmi$fitted, "data.frame")
  expect_is(geAmmi$trait, "character")
  expect_is(geAmmi$envMean, "numeric")
  expect_is(geAmmi$genoMean, "numeric")
  expect_is(geAmmi$overallMean, "numeric")
})

test_that("environmental scores are correct", {
  expect_equal(dim(geAmmi$envScores), c(3, 2))
  expect_equal(as.numeric(geAmmi$envScores),
               c(0.555925344086978, -0.795852871652761, 0.239927527565782,
                 -0.598008092310164, -0.182441424431847, 0.78044951674201))
})

test_that("genotypic scores are correct", {
  expect_equal(dim(geAmmi$genoScores), c(15, 2))
  expect_equal(as.numeric(geAmmi$genoScores),
               c(-12.9826230346155, 21.5755404989509, 5.86727127002252,
                 -20.9392111833447, 19.2865699416885, -16.283267546213,
                 6.47925849002732, -27.7204661464065, 2.63004569364296,
                 -0.446518670285875, 14.992858064254, -27.1548533454009,
                 -12.6486764594954, 1.37478052733848, 45.9692918998372,
                 -7.57810459752241, 9.5998190253813, -2.26762900778559,
                 -21.2510841390664, 2.80284705460814, 17.538092666485,
                 -1.39522497312836, 29.8807007034267, -10.2065925586386,
                 -10.4278856435621, 7.34435451407523, -21.7848037178221,
                 16.4865404970076, -5.88748313702092, -2.85354668643738))
})

test_that("importance is correct", {
  expect_equal(dim(geAmmi$importance), c(3, 3))
  expect_equal(as.numeric(as.matrix(geAmmi$importance)),
               c(20.3518912945909, 0.66873, 0.66873, 14.3243604621698, 0.33127,
                 1, 1.53416656938293e-14, 0, 1))
})

test_that("anova is correct", {
  expect_equal(dim(geAmmi$anova), c(6, 5))
  expect_equal(as.numeric(as.matrix(geAmmi$anova)),
               c(2, 14, 28, 15, 13, 0, 513.074183597498, 7644.67782855106,
                 8671.41494683825, 5798.79270973584, 2872.62223710241, 0,
                 256.537091798749, 546.048416325076, 309.693390958509,
                 386.586180649056, 220.97094131557, NaN, 0.828358303045345,
                 1.76319040789034, NA, NaN, NaN, NA, 0.447187338385732,
                 0.0980192661103099, NA, NaN, NaN, NA))
})

test_that("fitted values are correct", {
  expect_equal(dim(geAmmi$fitted), c(45, 3))
  expect_equal(geAmmi$fitted$fittedValue,
               c(69.9202841139079, 73.011566816012, 106.61526288738,
                 80.5465563933279, 82.5251967659839, 75.021141783374,
                 66.9476240572365, 62.7262039251742, 90.0880562676267,
                 97.9178080741426, 94.6684033674548, 68.598458437976,
                 49.8026925088618, 73.1067643214725, 126.941193543844,
                 83.1914449554764, 46.7063274979725, 96.6124038846183,
                 98.8912252179138, 56.4895429901647, 103.191496380399,
                 56.4800360998873, 111.48629093387, 81.162039649215,
                 93.0586621122849, 76.3241948942898, 95.1235008671379,
                 72.6228960205703, 67.6724657203917, 62.4858087891847,
                 70.1078833003919, 85.9579013086156, 108.166602116302,
                 64.400797374313, 86.825503225137, 110.873350657662,
                 69.5081486406389, 119.206304548313, 81.7188269956098,
                 90.2157232867707, 106.585788918718, 53.6811189560119,
                 83.0568584464682, 71.0879053477207, 115.01271913324))
})

test_that("means are correct", {
  expect_length(geAmmi$envMean, 3)
  expect_length(geAmmi$genoMean, 15)
  expect_equivalent(geAmmi$envMean,
                    c(81.2291475509183, 80.0998890675584, 87.7603621503942))
  expect_equivalent(geAmmi$genoMean,
                    c(74.406537456592, 68.5585985408667, 103.798089629433,
                      81.2795263285182, 75.2800809937618, 96.3619962738117,
                      64.3119362659209, 97.8062664691192, 84.3229743041505,
                      93.7307311577327, 92.5261290601542, 72.4676927537086,
                      68.4941489919668, 70.6223784631949, 101.479907155423))
  expect_equal(geAmmi$overallMean, 83.0297995896236)
})

test_that("options nPC functions properly", {
  geAmmi1 <- gxeAmmi(BLUEs, trait = "t1", nPC = 1)
  expect_equal(ncol(geAmmi1$envScores), 1)
  expect_equal(ncol(geAmmi1$genoScores), 1)
  expect_equal(geAmmi1$envScores, geAmmi$envScores[, 1, drop = FALSE])
  expect_equal(geAmmi1$genoScores, geAmmi$genoScores[, 1, drop = FALSE])
  expect_equal(geAmmi1$importance, geAmmi$importance)
  ## Third PC is very close to zero.
  expect_error(gxeAmmi(BLUEs, trait = "t1", nPC = 3),
               "should be smaller than")
  expect_warning(expect_error(gxeAmmi(BLUEsYear, trait = "t1", nPC = 3,
                                      byYear = TRUE),
                              "All years were skipped"),
                 "is larger than the number of trials")
})

test_that("making algorithm decide nPC functions properly", {
  geAmmi1 <- gxeAmmi(BLUEs, trait = "t1", nPC = NULL)
  expect_equal(ncol(geAmmi1$envScores), 2)
  expect_equal(ncol(geAmmi1$genoScores), 2)
  expect_equal(geAmmi1$envScores, geAmmi$envScores)
  expect_equal(geAmmi1$genoScore, geAmmi$genoScores)
  expect_equal(geAmmi1$importance, geAmmi$importance)
  ## Use year data but ignore year to use 6 environments.
  geAmmiYear1 <- gxeAmmi(BLUEsYear, trait = "t1", nPC = NULL)
  expect_equal(ncol(geAmmiYear1$envScores), 2)
})

test_that("option center functions properly", {
  geAmmiNC <- gxeAmmi(BLUEs, trait = "t1", center = FALSE)
  expect_equal(geAmmiNC$envScores, geAmmi$envScores)
  expect_equal(geAmmiNC$importance, geAmmi$importance)
  expect_equal(geAmmiNC$genoScores, geAmmi$genoScores)
})

test_that("option GGE functions properly", {
  geAmmiGGE <- gxeGGE(BLUEs, trait = "t1")
  expect_equal(dim(geAmmiGGE$anova), c(5, 5))
  expect_equal(as.numeric(as.matrix(geAmmiGGE$anova)),
               c(2, 42, 16, 14, 12, 513.074183597497, 16316.0927753893,
                 8142.70251527196, 5494.4162597362, 2678.97400038115,
                 256.537091798748, 388.478399414031, 508.918907204497,
                 392.458304266872, 223.247833365096, 0.660363850823368,
                 NA, 2.27961409315099, 1.75794899485117, NA, 0.521944696647847,
                 NA, 0.0769754584486909, 0.166762326653175, NA))
})

test_that("option excludeGeno functions properly", {
  expect_error(gxeAmmi(BLUEs, trait = "t1", excludeGeno = 1:10),
               "should be NULL or a character vector")
  expect_error(gxeAmmi(BLUEs, trait = "t1", excludeGeno = "g1"),
               "All genotypes to exclude should be in TD")
  expect_warning(gxeAmmi(BLUEs, trait = "t1", excludeGeno = paste0("G", 1:6)),
                 "Less than 10 genotypes present")
  geAmmi1 <- gxeAmmi(BLUEs, trait = "t1", excludeGeno = "G1")
  expect_equal(nrow(geAmmi1$genoScores), 14)
  expect_equal(nrow(geAmmi1$fitted), 42)
  expect_length(geAmmi1$genoMean, 14)
})

test_that("missing data is imputed correctly", {
  ## Add missing values.
  BLUEsMiss1 <- BLUEsMiss2 <- BLUEs
  BLUEsMiss1$E1[["t1"]] <- NA
  BLUEsMiss2$E1[1:5, "t1"] <- NA
  expect_error(gxeAmmi(BLUEsMiss1, trait = "t1"),
               "More than 30% missing values")
  geAmmiMiss2 <- gxeAmmi(BLUEsMiss2, trait = "t1")
  expect_equal(sum(is.na(geAmmiMiss2$genoScores)), 0)
  ## Add missing values for data by year.
  BLUEsYearMiss <- BLUEsYear
  BLUEsYearMiss$E1[["t1"]] <- BLUEsYearMiss$E2[["t1"]] <- NA
  expect_warning(geAmYear <- gxeAmmi(BLUEsYearMiss, trait = "t1",
                                     byYear = TRUE),
                 "More than 30% missing values")
  expect_is(geAmYear, "AMMI")
})

test_that("analysis is only run for years with at least three trials", {
  ## Delete E6 to leave only two trials for year 2.
  BLUEsYear[["E6"]] <- NULL
  expect_warning(gxeAmmi(BLUEsYear, trait = "t1", byYear = TRUE),
                 "less than 3 trials for 2")
  ## Delete E3 to leave only two trials for both years.
  BLUEsYear[["E3"]] <- NULL
  expect_warning(expect_error(gxeAmmi(BLUEsYear, trait = "t1", byYear = TRUE),
                              "All years were skipped"),
                 "less than 3 trials for 1")
})

geAmmiYear <- gxeAmmi(BLUEsYear, trait = "t1", byYear = TRUE)
test_that("output is of the right class when byYear = TRUE", {
  expect_is(geAmmiYear, "AMMI")
  expect_is(geAmmiYear$envScores, "list")
  expect_is(geAmmiYear$genoScores, "list")
  expect_is(geAmmiYear$importance, "list")
  expect_is(geAmmiYear$anova, "list")
  expect_is(geAmmiYear$fitted, "data.frame")
  expect_is(geAmmiYear$trait, "character")
  expect_is(geAmmiYear$envMean, "list")
  expect_is(geAmmiYear$genoMean, "list")
  expect_is(geAmmiYear$overallMean, "list")
})

test_that("output elements are of the right class when byYear = TRUE", {
  ## Only check two elements since function works equivalent for all elements.
  expect_length(geAmmiYear$envScores, 2)
  expect_named(geAmmiYear$envScores, c("1", "2"))
  expect_is(geAmmiYear$envScores[["1"]], "matrix")
  ## Different field names so not equal but equivalent.
  expect_length(geAmmiYear$anova, 2)
  expect_named(geAmmiYear$anova, c("1", "2"))
  expect_is(geAmmiYear$anova[["1"]], "data.frame")
})
