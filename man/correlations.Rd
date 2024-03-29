% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createVarComp.R
\name{correlations}
\alias{correlations}
\title{Compute different types of correlations.}
\usage{
correlations(varComp)
}
\arguments{
\item{varComp}{An object of class varComp.}
}
\value{
A list with three correlations.
}
\description{
Compute three types of correlations for models fitted with a nesting factor.
\itemize{
\item correlation between scenarios or environment types:
\deqn{\sigma_G^2 / (\sigma_G^2 + \sigma_{GS}^2)}
\item correlation between trials within scenarios or environment types:
\deqn{(\sigma_G^2 + \sigma_{GS}^2) / (\sigma_G^2 + \sigma_{GS}^2 +
\sigma_E^2)}
\item correlation between trials that belong to different
scenarios/environment types:
\deqn{\sigma_G^2 / (\sigma_G^2 + \sigma_{GS}^2 + \sigma_E^2)}
}
In these formulas the \eqn{\sigma} terms stand for the standard deviations of
the respective model terms. So \eqn{\sigma_S} is the standard deviation for
the scenario term in the model, \eqn{\sigma_{GS}} for the standard deviation
of the genotype by scenario term and \eqn{\sigma_E} corresponds to the
residual standard deviation.
}
\seealso{
Other Mixed model analysis: 
\code{\link{CRDR}()},
\code{\link{diagnostics}()},
\code{\link{gxeVarComp}()},
\code{\link{herit}()},
\code{\link{plot.varComp}()},
\code{\link{predict.varComp}()},
\code{\link{vc}()}
}
\concept{Mixed model analysis}
