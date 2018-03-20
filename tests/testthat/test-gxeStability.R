context("gxeStability")

testTD <- createTD(data = testData, genotype = "seed", trial = "field")

geStab <- gxeStability(TD = testTD, trait = "t1")
test_that("output is of the right class", {
  expect_is(geStab, "stability")
  expect_is(geStab$superiority, "data.frame")
  expect_is(geStab$static, "data.frame")
  expect_is(geStab$wricke, "data.frame")
  expect_identical(geStab$trait, "t1")
})

sup <- geStab$superiority
test_that("superiority is computed correctly", {
  expect_identical(dim(sup), c(15L, 3L))
  expect_equal(as.numeric(sup$genotype),
               c(7, 2, 14, 4, 12, 5, 9, 8, 10, 15, 3, 11, 13, 1, 6))
  expect_equal(sup$superiority,
               c(3940.5397803278, 3022.21904080722, 2598.17364359058, 2588.49105583387,
                 2518.79370648346, 2473.14071610774, 2429.67282734051, 2362.07050030342,
                 2189.94273125631, 2024.78044718482, 1986.14147273499, 1501.97635189109,
                 1422.63335221524, 1290.79686509336, 1179.03818107474))
  expect_equal(sup$mean,
               c(54.2031334320118, 66.2676045412444, 71.7251608269665, 73.0871656828908,
                 73.2435501712991, 75.3836554703585, 74.4525075789223, 74.8374571665672,
                 80.6390175697029, 81.4068278048538, 81.921276388438, 92.20470759419,
                 94.1532040202908, 93.6809432528475, 97.3247024025305))
})

sta <- geStab$static
test_that("static is computed correctly", {
  expect_identical(dim(sta), c(15L, 3L))
  expect_equal(as.numeric(sta$genotype),
               c(15, 6, 8, 13, 10, 7, 11, 5, 9, 4, 14, 3, 12, 2, 1))
  expect_equal(sta$static,
               c(1236.22441699745, 938.239475056893, 571.795429941847, 340.869668436155,
                 311.134712234077, 238.270628556254, 199.929473556735, 141.550914632088,
                 120.073390965768, 98.2382041653228, 85.1891664482231, 72.908961978567,
                 58.8949470780672, 53.127398526111, 3.53601858128704))
  expect_equal(sta$mean,
               c(81.4068278048538, 97.3247024025305, 74.8374571665672, 94.1532040202908,
                 80.6390175697029, 54.2031334320118, 92.20470759419, 75.3836554703585,
                 74.4525075789223, 73.0871656828908, 71.7251608269665, 81.921276388438,
                 73.2435501712991, 66.2676045412444, 93.6809432528475))
})

wri <- geStab$wricke
test_that("wricke is computed correctly", {
  expect_identical(dim(wri), c(15L, 3L))
  expect_equal(as.numeric(wri$genotype),
               c(15, 6, 10, 8, 13, 11, 5, 4, 9, 7, 3, 12, 14, 2, 1))
  expect_equal(wri$wricke,
               c(1902.30490526853, 1576.34662172865, 834.192141954125, 788.464006957686,
                 732.907499302012, 518.182498059668, 472.642532486111, 279.929271802651,
                 271.933243518989, 259.20458428707, 257.053039154622, 182.110068850864,
                 122.341672845375, 114.89546943427, 55.2006334326454))
  expect_equal(wri$mean,
               c(81.4068278048538, 97.3247024025305, 80.6390175697029, 74.8374571665672,
                 94.1532040202908, 92.20470759419, 75.3836554703585, 73.0871656828908,
                 74.4525075789223, 54.2031334320118, 81.921276388438, 73.2435501712991,
                 71.7251608269665, 66.2676045412444, 93.6809432528475))
})

geStabMin <- gxeStability(TD = testTD, trait = "t1", bestMethod = "min")
supMin <- geStabMin$superiority
test_that("option bestMethod for superiority is computed correctly", {
  expect_identical(dim(supMin), c(15L, 3L))
  expect_equal(as.numeric(supMin$genotype),
               c(6, 13, 1, 11, 15, 10, 3, 5, 9, 8, 12, 4, 14, 2, 7))
  expect_equal(supMin$superiority,
               c(1990.71087049962, 1670.54267372833, 1563.5969077264, 1549.85254273757,
                 1134.77672117742, 1112.14443245264, 1039.13064351809, 829.308210462886,
                 734.896555792387, 726.653499212872, 678.140843270401, 665.157598405587,
                 569.239996165416, 413.996898302837, 133.941993486244))
  expect_equal(supMin$mean,
               c(97.3247024025305, 94.1532040202908, 93.6809432528475, 92.20470759419,
                 81.4068278048538, 80.6390175697029, 81.921276388438, 75.3836554703585,
                 74.4525075789223, 74.8374571665672, 73.2435501712991, 73.0871656828908,
                 71.7251608269665, 66.2676045412444, 54.2031334320118))
})

staMin <- geStabMin$static
test_that("option bestMethod for static is computed correctly", {
  expect_identical(staMin, sta)
})

wriMin <- geStabMin$wricke
test_that("option bestMethod for wricke is computed correctly", {
  expect_identical(wriMin, wri)
})

geStabAsc <- gxeStability(TD = testTD, trait = "t1", sorted = "ascending")
test_that("option sort ascending functions properly", {
  expect_identical(geStabAsc$superiority$genotype, rev(sup$genotype))
  expect_identical(geStabAsc$superiority$superiority, rev(sup$superiority))
  expect_identical(geStabAsc$static$genotype, rev(sta$genotype))
  expect_identical(geStabAsc$static$static, rev(sta$static))
  expect_identical(geStabAsc$wricke$genotype, rev(wri$genotype))
  expect_identical(geStabAsc$wricke$wricke, rev(wri$wricke))
})

geStabNo <- gxeStability(TD = testTD, trait = "t1", sorted = "none")
test_that("option sort none functions properly", {
  expect_equal(as.numeric(geStabNo$superiority$genotype), 1:15)
  expect_identical(geStabNo$superiority$superiority,
                   sup$superiority[order(sup$genotype)])
  expect_equal(as.numeric(geStabNo$static$genotype), 1:15)
  expect_identical(geStabNo$static$static, sta$static[order(sta$genotype)])
  expect_equal(as.numeric(geStabNo$wricke$genotype), 1:15)
  expect_identical(geStabNo$wricke$wricke, wri$wricke[order(wri$genotype)])
})

