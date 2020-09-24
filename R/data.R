#' DROPS data set
#'
#' This dataset comes from the European Union project DROPS (DROught-tolerant
#' yielding PlantS). A panel of 256 maize hybrids was grown with two water
#' regimes (irrigated or rainfed), in seven fields in 2012 and 2013,
#' respectively, spread along a climatic transect from western to eastern
#' Europe, plus one site in Chile in 2013. This resulted in 28 experiments
#' defined as the combination of one year, one site and one water regime, with
#' two and three repetitions for rainfed and irrigated treatments, respectively.
#' A detailed environmental characterisation was carried out, with hourly
#' records of micrometeorological data and soil water status, and associated
#' with precise measurement of phenology. Grain yield and its components were
#' measured at the end of the experiment.\cr
#' 10 experiments have been selected from the full data set, two for each of
#' the five main environmental scenarios that were identified in the data. The
#' scenarios have been added to the data as well as a classification of the
#' genotypes in four genetic groups.\cr\cr
#'
#' The data.frame contains the genotypic means  (Best Linear Unbiased Estimators,
#' BLUEs), with one value per experiment  (Location × year × water regime) per
#' genotype.\cr
#' A data.frame with 2460 rows and 19 columns.\cr
#' \describe{
#' \item{Experiment}{experiments ID described by the three first letters of the
#' city’s name followed by the year of experiment and the water regime with W
#' for watered and R for rain-fed.}
#' \item{parent1}{identifier of donor dent line}
#' \item{Code_ID, Variety_ID, Accession_ID}{identifier of the genotype}
#' \item{geno.panel}{project in which the genetic material was generated}
#' \item{grain.yield}{genotypic mean for yield adjusted at 15\% grain moisture,
#' in ton per hectare (t ha^-1)}
#' \item{grain.number}{genotypic mean for number of grain per square meter}
#' \item{grain.weight}{genotypic mean for individual grain weight in milligram
#'  (mg)}
#' \item{anthesis}{genotypic mean for male flowering (pollen shed), in thermal
#' time cumulated since emergence (d_20°C)}
#' \item{silking}{genotypic mean for female flowering (silking emergence), in
#' thermal time cumulated since emergence (d_20°C)}
#' \item{plant.height}{genotypic mean for plant height, from ground level to
#' the base of the flag leaf (highest) leaf in centimeter (cm)}
#' \item{tassel.height}{genotypic mean for plant height including tassel, from
#' ground level to the highest point of the tassel in centimeter (cm)}
#' \item{ear.height}{genotypic mean for ear insertion height, from ground level
#' to ligule of the highest ear leaf in centimeter (cm)}
#' \item{year}{year in which the experiment was performed}
#' \item{loc}{location where the experiment was performed, a three letter
#' abbreviation}
#' \item{scenarioWater}{water scenario for the experiment, well watered (WW) or
#' water deficit (WD)}
#' \item{scenarioTemp}{temperature scenario for the experiment, Cool, Hot or
#' Hot(Day)}
#' \item{scenarioFull}{the full scenario for the experiment, a combination of
#' scenarioWater and scenarioTemp}
#' \item{geneticGroup}{the genetic group to which the genotype belongs}
#' }
#'
#' @source \url{https://data.inrae.fr/dataset.xhtml?persistentId=doi:10.15454/IASSTN}
#'
#' @references Millet, E. J., Pommier, C., et al. (2019). A multi-site
#' experiment in a network of European fields for assessing the maize yield
#' response to environmental scenarios (Data set).
#' \url{https://doi.org/10.15454/IASSTN}
"dropsPheno"

#' Field data for a maize experiment in Tlaltizapan, Mexico
#'
#' A dataset converted into a TD object containing data corresponding to an
#' F2 maize reference population from CIMMYT maize drought breeding program,
#' which was derived from the cross of a drought-tolerant line (P1) with a
#' drought susceptible line (P2) as described in detail by Ribaut et al. (1996,
#' 1997).\cr
#' DNA from 211 F2 plants was extracted to produce information for 132
#' co-dominant markers on 10 linkage groups. Phenotypic evaluations were
#' performed on 211 F2:3 families, each one derived from an original F2 plant.
#' The families were evaluated under different water and nitrogen regimes
#' during 1992, 1994 and 1996. In the winter of 1992 three water regimes were
#' imposed on the trials: well watered (WW), intermediate stress (IS) and severe
#' stress (SS). In the winter of 1994, only the IS and SS trials were available.
#' Nitrogen availability varied in the 1996 trials, with two low nitrogen
#' treatments (LN, in winter and summer) and one high-nitrogen treatment
#' (HN in summer). In each of the trials, five traits were evaluated but only
#' grain yield is included in the data.
#'
#' @format A TD object, a list containing 8 data.frames, each with the following
#' columns:
#' \describe{
#'   \item{trial}{trial, a combination of water regime, year and nitrogen
#'   treatment}
#'   \item{genotype}{genotype}
#'   \item{yld}{grain yield in tons}
#' }
#'
#' @source \url{https://link.springer.com/article/10.1007/BF00221905}
#'
#' @references Ribaut JM, Hoisington DA, Deutsch JA, Jiang C, Gonzalez de Leon D
#' (1996) Identification of quantitative trait loci under drought conditions in
#' tropical maize.1. Flowering parameters and the anthesis-silking interval.
#' Theor Appl Genet 92:905–914
#' @references Ribaut JM, Jiang C, Gonzalez de Leon D, Edmeades GO, Hoisington
#' DA (1997) Identification of quantitative trait loci under drought conditions
#' in tropical maize.2. Yield components and marker-assisted selection
#' strategies. Theor Appl Genet 94:887–896
"TDMaize"
