#' Field data for a wheat experiment in Chili.
#'
#' A dataset containing field data from an experiment with wheat in Chili
#' described in detail by Lado (2013). The experiment was performed in 2
#' locations in Chili with 2 different drought regimes in 2011 and 2012 for
#' the first location and 1 trial in 2012 for the second location. For 384
#' genotypes 4 traits were measured in 2011 but in 2012 only grain yield was
#' measured.
#'
#' @format A data frame with 4000 rows and 11 variables:
#' \describe{
#'   \item{rep}{replicate}
#'   \item{bl}{block id}
#'   \item{trt}{genotype}
#'   \item{row}{row within the field}
#'   \item{col}{column within the field}
#'   \item{DH}{Days to Heading, the number of days from sowing till 50% of
#'   the spikes emerged}
#'   \item{GY}{Grain Yield, in tons}
#'   \item{NKS}{Number of Kernels per Spike, calculated from 25 randomly
#'   selected spikes per plot}
#'   \item{TKW}{Thousand Kernel Weight, in grams, calculated from 25 randomly
#'   selected spikes per plot}
#'   \item{trial}{trail, a combination of location and year}
#'   \item{year}{year}
#' }
#'
#' @source \url{http://www.g3journal.org/content/3/12/2105/}
#' @references Lado, Bettina, Ivan Matus, Alejandra Rodríguez, Luis Inostroza,
#' Jesse Poland, François Belzile, Alejandro del Pozo, Martín Quincke,
#' Marina Castro, and Jarislav von Zitzewitz. 2013. “Increased Genomic
#' Prediction Accuracy in Wheat Breeding Through Spatial Adjustment of Field
#' Trial Data.” G3: Genes|Genomes|Genetics 3 (12): 2105–14.
#' doi:10.1534/g3.113.007807.
"wheatChl"

#' Test data
"TDHeat05"

#' Test data maize
"TDMaize"

#' Random test data for unit testing.
#' @keywords internal
"testData"

#' Fake data of f2 cross for unit testing.
#'
#' Base data taken from qtl package. Decreased size by removing most individuals
#' and SNPs. To be used for (unit) testing only.
#'
#' @keywords internal
"testF2"

#' Fake data of 4way cross for unit testing.
#'
#' Base data taken from qtl package. Decreased size by removing most individuals
#' and SNPs. To be used for (unit) testing only.
#'
#' @keywords internal
"test4way"

#' Fake data of backcross for unit testing.
#'
#' Base data taken from qtl package. Decreased size by removing most individuals
#' and SNPs. To be used for (unit) testing only.
#'
#' @keywords internal
"testBc"
