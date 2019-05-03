#' Field data for a wheat experiment in Chili.
#'
#' A dataset containing field data from an experiment with wheat in Chili
#' described in detail by Lado (2013). The experiment was performed in 2
#' locations in Chili with 2 different drought regimes in 2011 and 2012 for
#' the first location and 1 trial in 2012 for the second location. For 384
#' genotypes 4 traits were measured in 2011 but in 2012 only grain yield was
#' measured.
#'
#' @format A data.frame with 4000 rows and 11 columns:
#' \describe{
#'   \item{rep}{replicate}
#'   \item{bl}{block id}
#'   \item{trt}{genotype}
#'   \item{row}{row within the field}
#'   \item{col}{column within the field}
#'   \item{DH}{Days to Heading, the number of days from sowing till 50\% of
#'   the spikes emerged}
#'   \item{GY}{Grain Yield, in tons}
#'   \item{NKS}{Number of Kernels per Spike, calculated from 25 randomly
#'   selected spikes per plot}
#'   \item{TKW}{Thousand Kernel Weight, in grams, calculated from 25 randomly
#'   selected spikes per plot}
#'   \item{trial}{trial, a combination of location and year}
#'   \item{year}{year}
#' }
#'
#' @source \url{http://www.g3journal.org/content/3/12/2105/}
#'
#' @references Lado, Bettina, Ivan Matus, Alejandra Rodríguez, Luis Inostroza,
#' Jesse Poland, François Belzile, Alejandro del Pozo, Martín Quincke,
#' Marina Castro, and Jarislav von Zitzewitz. 2013. “Increased Genomic
#' Prediction Accuracy in Wheat Breeding Through Spatial Adjustment of Field
#' Trial Data.” G3: Genes|Genomes|Genetics 3 (12): 2105–14.
#' doi:10.1534/g3.113.007807.
"wheatChl"

#' Field data for a wheat experiment in Mexico.
#'
#' A dataset converted to a TD object containing raw plot data for one trial
#' from a series of wheat trials conducted in Mexico by CIMMYT. The different
#' trials took place under different regimes of irrigation and temperature,
#' there were 4 trials across two years, labelled as DRIP05, HEAT05, HEAT06 and
#' IRRI06. The TD object only contains the data for HEAT05. Within each trial,
#' a set of 167 progeny of a RIL (Recombinant Inbred Line; 8 generations)
#' population were tested alongside the population parents (Seri and Babax). A
#' lattice design with two replicates was used for each trial. In the first
#' replicate the entries were not randomized, as they were considered to be a
#' random selection from a population.
#'
#' @format A TD object, a list containing 1 data.frames with the following
#' columns:
#' \describe{
#'   \item{trial}{trial, a combination of watering regime, year and nitrogen
#'   treatment}
#'   \item{genotype}{genotype}
#'   \item{Plot}{plot number in the field}
#'   \item{repId}{replicate}
#'   \item{subBlock}{block id}
#'   \item{rowId}{row within the field (as factor)}
#'   \item{colId}{column within the field (as factor)}
#'   \item{yield}{yield in grams per square meter}
#'   \item{rowCoord}{row within the field (as numerical value)}
#'   \item{colCoord}{column within the field (as numerical value)}
#' }
"TDHeat05"

#' Field data for a maize experiment in Tlaltizapan, Mexico.
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
#'   \item{trial}{trial, a combination of watering regime, year and nitrogen
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

#' Random test data for unit testing.
#' @keywords internal
"testData"
