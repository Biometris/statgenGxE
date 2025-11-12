# S3 class AMMI

Function for creating objects of S3 class AMMI.  
[`print`](https://rdrr.io/r/base/print.html),
[`summary`](https://rdrr.io/r/base/summary.html),
[`plot`](https://rdrr.io/r/graphics/plot.default.html) and
[`report`](https://rdrr.io/pkg/statgenSTA/man/report.html) methods are
available.

## Usage

``` r
createAMMI(
  envScores,
  genoScores,
  importance,
  anova,
  fitted,
  trait,
  envMean,
  genoMean,
  overallMean,
  dat,
  GGE,
  byYear
)
```

## Arguments

- envScores:

  A matrix containing environmental scores.

- genoScores:

  A matrix containing genotypic scores.

- importance:

  A data.frame containing the importance of the principal components.

- anova:

  A data.frame containing anova scores of the AMMI analysis.

- fitted:

  A matrix containing fitted values from the AMMI model.

- trait:

  A character string indicating the analyzed trait.

- envMean:

  A numerical vector containing the environmental means.

- genoMean:

  A numerical vector containing the genotypic means.

- overallMean:

  A numerical value containing the overall mean.

- GGE:

  Has a GGE analysis been performed?

- byYear:

  Has the analysis been performed by year?

## See also

[`plot.AMMI`](plot.AMMI.md), [`report.AMMI`](report.AMMI.md)
