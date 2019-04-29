context("ExtractAsreml")

if (requireNamespace("asreml", quietly = TRUE)) {
  testTD <- createTD(data = testData[testData$field == "E1", ],
                     genotype = "seed", repId = "rep",
                     subBlock = "block", rowId = "Y", colId = "X",
                     rowCoord = "Y", colCoord = "X")

  modelAs <- fitTD(testTD, design = "res.ibd", traits = "t1", engine = "asreml")

  test_that("the output of extract is of the proper type - asreml", {
    expect_is(STExtract(modelAs, what = "BLUEs"), "list")
    expect_is(STExtract(modelAs), "list")
    expect_length(STExtract(modelAs)[[1]], 23)
    expect_is(STExtract(modelAs, what = c("BLUEs", "BLUPs")), "list")
    expect_length(STExtract(modelAs, what = c("BLUEs", "BLUPs"))[[1]], 2)
  })

  extAs <- STExtract(modelAs)[[1]]
  test_that("BLUEs are computed correctly", {
    expect_is(extAs$BLUEs, "data.frame")
    expect_identical(dim(extAs$BLUEs), c(15L, 2L))
    expect_named(extAs$BLUEs, c("genotype", "t1"))
    expect_equal(extAs$BLUEs$t1,
                 c(91.8559516394161, 73.6623287954585, 74.8843349808429,
                   84.5319987913184, 68.9075264154031, 65.0344131821085,
                   55.3031580266259, 77.0436561256374, 62.1348724815413,
                   67.3983425459108, 108.495613219885, 65.3049181280563,
                   114.94121644052, 78.9340989737003, 76.3592209421872))
  })

  test_that("SE of BLUEs are computed correctly", {
    expect_is(extAs$seBLUEs, "data.frame")
    expect_identical(dim(extAs$seBLUEs), c(15L, 2L))
    expect_named(extAs$seBLUEs, c("genotype", "t1"))
    expect_equal(extAs$seBLUEs$t1, rep(x = 14.5728200113826, times = 15))
  })

  test_that("BLUPs are computed correctly", {
    expect_is(extAs$BLUPs, "data.frame")
    expect_identical(dim(extAs$BLUPs), c(15L, 2L))
    expect_named(extAs$BLUPs, c("genotype", "t1"))
    expect_equal(extAs$BLUPs$t1,
                 c(83.1254685686293, 75.9302650186231, 76.6335759637273,
                   78.3870979895003, 75.3085056773248, 73.8785718154553,
                   68.9076125899214, 75.5226969514441, 72.7178660892812,
                   74.782809512082, 89.2545078476149, 72.9693083410864,
                   90.019070308295, 79.1438169239316, 78.2104770916953))
  })

  test_that("SE of BLUPs are computed correctly", {
    expect_is(extAs$seBLUPs, "data.frame")
    expect_identical(dim(extAs$seBLUPs), c(15L, 2L))
    expect_named(extAs$seBLUPs, c("genotype", "t1"))
    expect_equal(extAs$seBLUPs$t1, rep(x = 8.96361730231537, times = 15))
  })

  test_that("unit errors are computed correctly", {
    expect_is(extAs$ue, "data.frame")
    expect_identical(dim(extAs$ue), c(15L, 2L))
    expect_named(extAs$ue, c("genotype", "t1"))
    expect_equal(extAs$ue$t1, rep(x = 195.688651981688, times = 15))
  })

  test_that("heritability is computed correctly", {
    expect_is(extAs$heritability, "numeric")
    expect_length(extAs$heritability, 1)
    expect_named(extAs$heritability, "t1")
    expect_equivalent(extAs$heritability, 0.34526827189411)
  })

  test_that("varGen is computed correctly", {
    expect_is(extAs$varGen, "numeric")
    expect_length(extAs$varGen, 1)
    expect_named(extAs$varGen, "t1")
    expect_equivalent(extAs$varGen, 99.3045736226269)
  })

  test_that("varErr is computed correctly", {
    expect_is(extAs$varErr, "numeric")
    expect_length(extAs$varErr, 1)
    expect_named(extAs$varErr, "t1")
    expect_equivalent(extAs$varErr, 320.639866476272)
  })

  test_that("fitted values are computed correctly", {
    expect_is(extAs$fitted, "data.frame")
    expect_identical(dim(extAs$fitted), c(30L, 3L))
    expect_named(extAs$fitted, c("genotype", "repId", "t1"))
    expect_equal(extAs$fitted$t1,
                 c(68.7325142281377, 63.8332414120377, 65.1408096443281,
                   95.2835477394968, 48.3415831521691, 61.4693120482354,
                   75.9281618109135, 93.7839958782487, 105.68921935359,
                   124.19321352745, 82.7008157447754, 55.1142370860309,
                   70.963443679784, 88.4283555393336, 72.2011477081544,
                   86.2956532125677, 78.3119310809243, 53.8419769393217,
                   92.7273883030725, 68.5995143159817, 72.7941198083142,
                   79.9243220760605, 67.7916590387072, 71.4567388807611,
                   75.1235098827627, 56.7643391139301, 75.2800017043882,
                   109.956794307189, 107.034432132581, 61.8773220279744))
  })

  test_that("residuals are computed correctly", {
    expect_is(extAs$resid, "data.frame")
    expect_identical(dim(extAs$resid), c(30L, 3L))
    expect_named(extAs$resid, c("genotype", "repId", "t1"))
    expect_equal(extAs$resid$t1,
                 c(-3.73900390681584, 6.66109558076825, -25.5878773515572,
                   9.55404865955619, -10.0401315668725, -2.71192069907498,
                   10.0401315668725, -7.36190314933543, -10.4901112979641,
                   10.490111297964, -8.8681127063347, 8.86811270633465,
                   -6.66109558076828, -9.55404865955656, 20.7840323904322,
                   7.41427077513407, -16.0720951734068, 0.554803295555431,
                   25.5878773515572, 2.71192069907495, 5.81683105357236,
                   -5.81683105357239, -7.41427077513409, 16.0720951734068,
                   -20.7840323904322, -0.554803295555459, 7.3619031493354,
                   4.05951690606135, -4.05951690606138, 3.73900390681578))
  })

  test_that("standardized residuals are computed correctly", {
    expect_is(extAs$stdRes, "data.frame")
    expect_identical(dim(extAs$stdRes), c(30L, 3L))
    expect_named(extAs$stdRes, c("genotype", "repId", "t1"))
    expect_equal(extAs$stdRes$t1,
                 c(-0.322220969016617, 0.574041837408854, -2.20511955611155,
                   0.823351591449873, -0.865241385983749, -0.233708692831061,
                   0.865241385983746, -0.634436236416103, -0.904019870468839,
                   0.904019870468837, -0.764238802846611, 0.764238802846607,
                   -0.574041837408857, -0.823351591449903, 1.79113240419721,
                   0.638949188943159, -1.3850657046552, 0.0478120001911751,
                   2.20511955611155, 0.233708692831059, 0.501284562787259,
                   -0.501284562787262, -0.63894918894316, 1.3850657046552,
                   -1.79113240419721, -0.0478120001911775, 0.6344362364161,
                   0.349842231730744, -0.349842231730747, 0.322220969016612))
  })

  test_that("rMeans are computed correctly", {
    expect_is(extAs$rMeans, "data.frame")
    expect_identical(dim(extAs$rMeans), c(30L, 3L))
    expect_named(extAs$rMeans, c("genotype", "repId", "t1"))
    expect_equal(extAs$rMeans$t1,
                 c(76.2962530364337, 68.4670155927884, 62.7365586150466,
                   86.4476589145487, 56.3136159266827, 67.5632011921593,
                   83.8944436938505, 92.5961231482435, 85.7209217331892,
                   104.222650253522, 86.4838705395223, 58.9030427723545,
                   75.5993166445069, 79.5902624014691, 74.9806082418925,
                   89.7330630070549, 79.9588053249045, 67.9612432902301,
                   90.3173863822144, 74.6955022438778, 71.893078596555,
                   79.0253796482735, 71.2313344867217, 73.1014088118249,
                   77.9066838647156, 70.8873189130533, 74.0943946279102,
                   91.2246892865721, 88.2986136637489, 69.4388565233541))
  })

  test_that("random effects are computed correctly", {
    expect_is(extAs$ranEf, "data.frame")
    expect_identical(dim(extAs$ranEf), c(15L, 2L))
    expect_named(extAs$ranEf, c("genotype", "t1"))
    expect_equal(extAs$ranEf$t1,
                 c(5.47011701726975, -1.72164155135825, -1.01873657237448,
                   0.734594632206112, -2.343499215799, -3.77278047184219,
                   -8.74100650302057, -2.12846550898244, -4.93292606147077,
                   -2.86896607121313, 11.5963638704982, -4.68128886084527,
                   12.3611217374851, 1.49001662689314, 0.557096932553493))
  })

  test_that("residuals are computed correctly for genotype random", {
    expect_is(extAs$residR, "data.frame")
    expect_identical(dim(extAs$residR), c(30L, 3L))
    expect_named(extAs$residR, c("genotype", "repId", "t1"))
    expect_equal(extAs$residR$t1,
                 c(-11.3027427151119, 2.02732140001756, -23.1836263222757,
                   18.3899374845042, -18.0121643413861, -8.80580984299887,
                   2.07384968393552, -6.17403041933025, 9.47818632243666,
                   30.460674571892, -12.6511675010816, 5.07930702001109,
                   -11.2969685454911, -0.715955521692109, 18.0045718566941,
                   3.97686098064682, -17.718969417387, -13.564463055353,
                   27.9978792724153, -3.38406722882111, 6.71787226533151,
                   -4.91788862578539, -10.8539462231486, 14.427425242343,
                   -23.5672063723852, -14.6777830946787, 8.54751022581337,
                   22.7916219266781, 14.6763015627702, -3.82253058856389))
  })

  test_that("standardized residuals are computed correctly for genotype random", {
    expect_is(extAs$stdResR, "data.frame")
    expect_identical(dim(extAs$stdResR), c(30L, 3L))
    expect_named(extAs$stdResR, c("genotype", "repId", "t1"))
    expect_equal(extAs$stdResR$t1,
                 c(-0.77349533987043, 0.1387383305856, -1.58655534975611,
                   1.25850258678844, -1.23264994437278, -0.602619474672998,
                   0.141922461348357, -0.422515479467228, 0.648633091596991,
                   2.08455508760423, -0.865773852652894, 0.347598844703141,
                   -0.773100188586724, -0.0489959183838563, 1.23213035796129,
                   0.272153716436463, -1.21258535358246, -0.928274598972526,
                   1.91601540345721, -0.231586288149483, 0.45973291811122,
                   -0.336552289114823, -0.74278226395305, 0.987330816294779,
                   -1.61280538381468, -1.0044638818668, 0.584942919944722,
                   1.55972880147079, 1.00436249426062, -0.261592223347099))
  })

  test_that("Waldtest is computed correctly", {
    expect_is(extAs$wald$t1$Ftest, "numeric")
    expect_length(extAs$wald$t1$Ftest, 4)
    expect_equivalent(unlist(extAs$wald$t1$Ftest),
                      c(1.46, 14, 14, 0.244132249447016))
  })

  test_that("CV is computed correctly", {
    expect_is(extAs$CV, "numeric")
    expect_length(extAs$CV, 1)
    expect_named(extAs$CV, "t1")
    expect_equivalent(extAs$CV, 23.4580700451519)
  })

  test_that("rDf is computed correctly", {
    expect_is(extAs$rDf, "integer")
    expect_length(extAs$rDf, 1)
    expect_named(extAs$rDf, "t1")
    expect_equivalent(extAs$rDf, 14)
  })

  test_that("rDfR is computed correctly", {
    expect_is(extAs$rDfR, "integer")
    expect_length(extAs$rDfR, 1)
    expect_named(extAs$rDfR, "t1")
    expect_equivalent(extAs$rDfR, 28)
  })

  test_that("sed is computed correctly", {
    expect_is(extAs$sed$t1, "numeric")
    expect_length(extAs$sed$t1, 3)
    expect_equivalent(extAs$sed$t1,
                      c(18.215842753241, 20.2671887094224, 20.6090797021193))
  })

  test_that("lsd is computed correctly", {
    expect_is(extAs$lsd$t1, "numeric")
    expect_length(extAs$lsd$t1, 3)
    expect_equivalent(extAs$lsd$t1,
                      c(39.0690970463554, 43.4687965454872, 44.2020797953426))
  })

  test_that("correct attributes are added", {
    expect_equal(attr(x = extAs, which = "traits"), "t1")
    expect_equal(attr(x = extAs, which = "design"), "res.ibd")
    expect_equal(attr(x = extAs, which = "engine"), "asreml")
  })

  test_that("calculated values are logically correct", {
    ## Fitted + residuals should match raw data.
    expect_equal(testTD[[1]]$t1, extAs$fitted$t1 + extAs$resid$t1)
    expect_equal(testTD[[1]]$t1, extAs$rMeans$t1 + extAs$residR$t1)
  })
}

