# DROPS data set

This dataset comes from the European Union project DROPS
(DROught-tolerant yielding PlantS). A panel of 256 maize hybrids was
grown with two water regimes (irrigated or rainfed), in seven fields in
2012 and 2013, respectively, spread along a climatic transect from
western to eastern Europe, plus one site in Chile in 2013. This resulted
in 28 experiments defined as the combination of one year, one site and
one water regime, with two and three repetitions for rainfed and
irrigated treatments, respectively. A detailed environmental
characterisation was carried out, with hourly records of
micrometeorological data and soil water status, and associated with
precise measurement of phenology. Grain yield and its components were
measured at the end of the experiment.  
10 experiments have been selected from the full data set, two for each
of the five main environmental scenarios that were identified in the
data. The scenarios have been added to the data as well as a
classification of the genotypes in four genetic groups.  
  

## Usage

``` r
dropsPheno
```

## Format

An object of class `data.frame` with 2460 rows and 20 columns.

## Source

[doi:10.15454/IASSTN](https://doi.org/10.15454/IASSTN)

## Details

The data.frame contains the genotypic means (Best Linear Unbiased
Estimators, BLUEs), with one value per experiment (Location × year ×
water regime) per genotype.  
A data.frame with 2460 rows and 19 columns.  

- Experiment:

  experiments ID described by the three first letters of the city’s name
  followed by the year of experiment and the water regime with W for
  watered and R for rain-fed.

- parent1:

  identifier of donor dent line

- Code_ID, Variety_ID, Accession_ID:

  identifier of the genotype

- geno.panel:

  project in which the genetic material was generated

- grain.yield:

  genotypic mean for yield adjusted at 15\\ in ton per hectare (t ha^-1)

- grain.number:

  genotypic mean for number of grain per square meter

- grain.weight:

  genotypic mean for individual grain weight in milligram (mg)

- anthesis:

  genotypic mean for male flowering (pollen shed), in thermal time
  cumulated since emergence (d_20°C)

- silking:

  genotypic mean for female flowering (silking emergence), in thermal
  time cumulated since emergence (d_20°C)

- plant.height:

  genotypic mean for plant height, from ground level to the base of the
  flag leaf (highest) leaf in centimeter (cm)

- tassel.height:

  genotypic mean for plant height including tassel, from ground level to
  the highest point of the tassel in centimeter (cm)

- ear.height:

  genotypic mean for ear insertion height, from ground level to ligule
  of the highest ear leaf in centimeter (cm)

- year:

  year in which the experiment was performed

- loc:

  location where the experiment was performed, a three letter
  abbreviation

- scenarioWater:

  water scenario for the experiment, well watered (WW) or water deficit
  (WD)

- scenarioTemp:

  temperature scenario for the experiment, Cool, Hot or Hot(Day)

- scenarioFull:

  the full scenario for the experiment, a combination of scenarioWater
  and scenarioTemp

- geneticGroup:

  the genetic group to which the genotype belongs

## References

Millet, E. J., Pommier, C., et al. (2019). A multi-site experiment in a
network of European fields for assessing the maize yield response to
environmental scenarios (Data set).
[doi:10.15454/IASSTN](https://doi.org/10.15454/IASSTN)
