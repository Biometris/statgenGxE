context("QTLDetect")

QTLDetF2 <- QTLDetect(testF2, trait = "phenotype")
QTLDet4way <- QTLDetect(test4way, trait = "phenotype")
QTLDetBc <- QTLDetect(testBc, trait = "pheno1")
test_that("QTLDetect on different types of cross objects functions properly", {
  expect_is(QTLDetF2, "QTLDet")
  expect_is(QTLDet4way, "QTLDet")
  expect_is(QTLDetBc, "QTLDet")
  expect_error(QTLDetect(test4way, trait = "phenotype", type = "CIM"),
               "CIM has not been implemented")
})

test_that("QTLDetect gives correct scores for different types", {
  expect_equal(QTLDetF2$scores$lod,
               c(1.40144553960757, 4.3656024143265, 2.48961022964713,
                 1.28380587618962, 0.971974741439599, 0.528546101782021,
                 0.279985665897152, 0.279985665897152, 0.199644978888341,
                 1.84895395868088, 1.83759097042948, 0.394673215857738,
                 0.0811785321388383, 1.22278473130343, 2.42361178806114,
                 0.143653670016825))
  expect_equal(QTLDetF2$peaks$altName, c("Q1@9.8", "QX@14.2"))
  expect_equal(rownames(QTLDetF2$peaks), c("D1M318", "DXM66"))
  QTLDetF2S <- QTLDetect(testF2, trait = "phenotype", type = "SIM")
  expect_equal(QTLDetF2S$peaks$altName, "Q1@9.8")
  expect_equal(rownames(QTLDetF2S$peaks), "D1M318")
  #QTLDetF2C <- QTLDetect(testF2, trait = "phenotype", type = "CIM")
  #expect_equal(rownames(QTLDetF2C$peaks), c("D1M318", "c2.loc30"))
})

test_that("option thrType functions properly", {
  QTLDetThrL <- QTLDetect(testF2, trait = "phenotype", thrType = "liji")
  QTLDetThrB <- QTLDetect(testF2, trait = "phenotype", thrType = "bonferroni")
  QTLDetThrF <- QTLDetect(testF2, trait = "phenotype", thrType = "fixed")
  expect_equal(QTLDetThrL$info$thrType, "liji")
  expect_equal(QTLDetThrL$info$thr, 2.30102999566398)
  expect_equal(QTLDetThrB$info$thrType, "bonferroni")
  expect_equal(QTLDetThrB$info$thr, 3.07188200730613)
  expect_equal(QTLDetThrF$info$thrType, "fixed")
  expect_equal(QTLDetThrF$info$thr, 3)
})

test_that("option thrAlpha functions properly", {
  QTLDetThrL <- QTLDetect(testF2, trait = "phenotype", thrType = "liji",
                           thrAlpha = 0.2)
  QTLDetThrB <- QTLDetect(testF2, trait = "phenotype", thrType = "bonferroni",
                          thrAlpha = 0.2)
  expect_equal(QTLDetThrL$info$thr, 1.69897000433602)
  expect_equal(QTLDetThrB$info$thr, 2.46982201597816)
})

test_that("option thrDist functions properly", {
  QTLDetThrB <- QTLDetect(testF2, trait = "phenotype", thrType = "bonferroni",
                          thrDist = 8)
  expect_equal(QTLDetThrB$info$thr, 2.77815125038364)
})

test_that("option thrFixed functions properly", {
  QTLDetThr <- QTLDetect(testF2, trait = "phenotype", thrType = "fixed",
                         thrFixed = 2)
  expect_equal(rownames(QTLDetThr$peaks), c("D1M318", "DXM66"))
  QTLDetSThr <- QTLDetect(testF2, trait = "phenotype", type = "SIM",
                          thrType = "fixed", thrFixed = 2)
  expect_equal(rownames(QTLDetSThr$peaks), c("D1M318", "D2M241"))
  #QTLDetCThr <- QTLDetect(testF2, trait = "phenotype", type = "CIM",
  #                        thrType = "fixed", thrFixed = 2.5)
  #expect_equal(rownames(QTLDetCThr$peaks), c("c2.loc30"))
})

test_that("option window functions properly", {
  QTLDetWin <- QTLDetect(testF2, trait = "phenotype", thrType = "fixed",
                         thrFixed = 2, window = 2)
  expect_equal(rownames(QTLDetWin$peaks), c("D1M318", "D1M212", "DXM66"))
  QTLDetSWin <- QTLDetect(testF2, trait = "phenotype", type = "SIM",
                          thrType = "fixed", thrFixed = 2, window = 2)
  expect_equal(rownames(QTLDetSWin$peaks),
               c("c1.loc5", "D1M318", "c1.loc15", "c1.loc20", "c1.loc25",
                 "c1.loc30", "D2M241", "c2.loc25"))
})

test_that("option step functions properly", {
  QTLDetSt1 <- QTLDetect(testF2, trait = "phenotype", step = 1)
  QTLDetSt2 <- QTLDetect(testF2, trait = "phenotype", step = 50)
  expect_equal(dim(QTLDetSt1$cross$geno$`1`$prob), c(50, 145, 3))
  expect_equal(dim(QTLDetSt2$cross$geno$`1`$prob), c(50, 10, 3))
})
