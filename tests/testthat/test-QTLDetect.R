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
  expect_equal(QTLDetF2$peaks$altName, "Q1@9.8")
  expect_equal(rownames(QTLDetF2$peaks), "D1M318")
  QTLDetF2S <- QTLDetect(testF2, trait = "phenotype", type = "SIM")
  expect_equal(nrow(QTLDetF2S$peaks), 0)
  ## CIM sometimes seems to give random warnings.
  #QTLDetF2C <- QTLDetect(testF2, trait = "phenotype", type = "CIM")
  #expect_equal(nrow(QTLDetF2C$peaks), 0)
})

test_that("option thr functions properly", {
  QTLDetThr <- QTLDetect(testF2, trait = "phenotype", thr = 2)
  expect_equal(rownames(QTLDetThr$peaks), c("D1M318", "DXM66"))
  QTLDetSThr <- QTLDetect(testF2, trait = "phenotype", type = "SIM", thr = 2)
  expect_equal(rownames(QTLDetSThr$peaks), c("D1M318", "D2M241"))
  ## CIM sometimes seems to give random warnings.
  #QTLDetCThr <- QTLDetect(testF2, trait = "phenotype", type = "CIM", thr = 2.5)
  #expect_equal(rownames(QTLDetCThr$peaks), c("c1.loc20", "c2.loc20"))
})

test_that("option window functions properly", {
  QTLDetWin <- QTLDetect(testF2, trait = "phenotype", thr = 2, window = 2)
  expect_equal(rownames(QTLDetWin$peaks), c("D1M318", "D1M212", "DXM66"))
  QTLDetSWin <- QTLDetect(testF2, trait = "phenotype", type = "SIM", thr = 2,
                          window = 2)
  expect_equal(rownames(QTLDetSWin$peaks),
               c("c1.loc5", "D1M318", "c1.loc15", "c1.loc20", "c1.loc25",
                 "c1.loc30", "D2M241", "c2.loc25"))
})
