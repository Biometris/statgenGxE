% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gxeVarCov.R
\name{gxeVarCov}
\alias{gxeVarCov}
\title{Selects the best variance-covariance model for a set of trials}
\usage{
gxeVarCov(
  TD,
  trials = names(TD),
  trait,
  models = c("identity", "cs", "diagonal", "hcs", "outside", "fa", "fa2", "unstructured"),
  engine = c("lme4", "asreml"),
  criterion = c("BIC", "AIC"),
  ...
)
}
\arguments{
\item{TD}{An object of class \code{\link[statgenSTA]{TD}}.}

\item{trials}{A character string specifying the trials to be analyzed. If
not supplied, all trials are used in the analysis.}

\item{trait}{A character string specifying the trait to be analyzed.}

\item{models}{A character vector specifying the models to be fitted.}

\item{engine}{A character string specifying the engine used for modeling.
Either "lme4" or "asreml".}

\item{criterion}{A string specifying a goodness-of-fit criterion. Either
"AIC" or "BIC".}

\item{...}{Further arguments to be passed to the modeling engine.}
}
\value{
An object of class \code{\link{varCov}}, a list object containing:
\item{STA}{An object of class STA containing the best fitted model.}
\item{choice}{A character string indicating the best fitted model.}
\item{summary}{A data.frame with a summary of the fitted models.}
\item{vcov}{The covariance matrix of the best fitted model.}
\item{criterion}{A character string indicating the goodness-of-fit criterion
used for determining the best model, either "AIC" or "BIC".}
\item{engine}{A character string containing the engine used for
the analysis.}
\item{trait}{A character string containing the trait analyzed.}
\item{dat}{A data.frame with the full data set used for the analysis.}
}
\description{
This function selects the best covariance structure for genetic correlations
between trials. It fits a range of variance-covariance models (identity,
compound symmetry (cs), diagonal, simple correlation with heterogeneous
variance (outside), heterogeneous compound symmetry (hcs),
first order factor analytic (fa), second order factor analytic (fa2) and
unstructured), and selects the best one using a goodness-of-fit criterion.
See details for the exact models fitted.
}
\details{
The models fitted are of the form \eqn{y_{ij} = \mu_j + \epsilon_{ij}}, where
\eqn{y_{ij}} is the phenotypic value of genotype \eqn{i} in environment
\eqn{j}, \eqn{\mu_j}
is the environmental mean, and \eqn{\epsilon_{ij}} represents mainly genetic
variation, although some non-genetic variation may be included as well.
The random term \eqn{\epsilon_{ij}} is modeled in eight ways as described in
the table below.\tabular{lllll}{
   Model \tab Description \tab var(\eqn{g_{ij}}) \tab cov(\eqn{g_{ij}};\eqn{g_{ik}}) \tab Number of parameters \cr
   identity \tab identity \tab \eqn{\sigma_G^2} \tab 0 \tab 1 \cr
   cs \tab compound symmetry \tab \eqn{\sigma_G^2+\sigma_{GE}^2} \tab \eqn{\sigma_{GE}^2} \tab 2 \cr
   diagonal \tab diagonal matrix (heteroscedastic) \tab \eqn{\sigma_{GE_j}^2} \tab 0 \tab \eqn{J} \cr
   hcs \tab heterogeneous compound symmetry \tab \eqn{\sigma_G^2+\sigma_{GE_j}^2} \tab \eqn{\sigma_G^2} \tab \eqn{J+1} \cr
   outside \tab heterogeneity outside \tab \eqn{\sigma_{G_j}^2} \tab \eqn{\theta} \tab \eqn{J+1} \cr
   fa \tab first order factor analytic \tab \eqn{\lambda_{1j}^2+\sigma_{GE_j}^2} \tab \eqn{\lambda_{1j}\lambda_{1k}} \tab \eqn{2J} \cr
   fa2 \tab second order factor analytic \tab \eqn{\lambda_{1j}^2+\lambda_{2j}^2+\sigma_{GE_j}^2} \tab \eqn{\lambda_{1j}\lambda_{1k}+\lambda_{2j}\lambda_{2k}} \tab \eqn{3J-1} \cr
   unstructured \tab unstructured \tab \eqn{\sigma_{G_j}^2} \tab \eqn{\sigma_{G_{j,k}}^2} \tab \eqn{J(J+1)/2} \cr
}


In this table \eqn{J} is the number of environments, \eqn{\sigma_G^2} the
variance component for the genotype main effects, \eqn{\sigma_{GE}^2} the
variance component for GxE interactions. \eqn{\sigma_{G_j}^2} and
\eqn{\sigma_{GE_j}^2} are the environment specific variance components for
the genotype main effects and GxE interaction in environment \eqn{j}.
\eqn{\sigma_{G_{j,k}}^2} is the genetic covariance between environments
\eqn{j} and \eqn{k}. \eqn{\theta} is the common correlation between
environments and \eqn{\lambda_{1j}} and \eqn{\lambda_{2j}} are
environment specific multiplicative parameters.
}
\note{
If \code{engine = "lme4"}, only the compound symmetry model can be
fitted.
}
\examples{
## Select the best variance-covariance model using lme4 for modeling.
geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld")

## Summarize results.
summary(geVarCov)

\donttest{
## Create a pdf report summarizing the results.
report(geVarCov, outfile = tempfile(fileext = ".pdf"))
}

\donttest{
if (requireNamespace("asreml", quietly = TRUE)) {
  ## Select the best variance-covariance model using asreml for modeling.
  ## Use BIC as a goodness-of-fit criterion.
  geVarCov2 <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml",
                         criterion = "BIC")

  summary(geVarCov2)

  ## Plot a heatmap of the correlation matrix for the best model.
  plot(geVarCov2)
  }
}

}
\seealso{
Other varCov: 
\code{\link{fitted.varCov}()},
\code{\link{plot.varCov}()},
\code{\link{report.varCov}()},
\code{\link{residuals.varCov}()}
}
\concept{varCov}
