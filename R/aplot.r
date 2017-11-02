#' Residual diagnostic plots
#'
#' This function is to produce a histogram of residuals, a normal Q-Q plot, a residuals vs fitted values plot,
#' and an absolute residuals vs fitted values plot.
#' @param x A fitted model object.
#' @param ... Other graphical parameters (see \code{\link[lattice]{xyplot}} for details).
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Row","Column"),
#'                      trait.names="yield", env ="Env", rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.run.model(data=mydat, design="res.rowcol", trait="yield",
#'                         genotype="Genotype", rep="Rep", row="Row", col="Column",
#'                         tryspatial="always")
#' aplot(mymodel$mfix)
#' #aplot(mymodel$mmix)
#' #c.f. in-built plot for "asreml" or "lme4"
#' #plot(mymodel$mfix)
#'
#' @export
aplot <- function(x, ...) {

  # Diagnostic plots
  if (class(x)=="asreml"){
    Resid <- x$residuals
    Fitted <- x$fitted.values
  }else{
    Resid <- residuals(x)
    Fitted <- fitted(x)
  }

  trellis.obj <- vector("list", 4)
  names(trellis.obj) <- c("histogram","qq", "ResidFitted", "AbsResidFitted")

  # Histogram of residuals
  trellis.obj[["histogram"]] <- histogram(x=~Resid, xlab="Residuals", ...)

  # Q-Q plot of residuals
  trellis.obj[["qq"]] <- qqmath(~Resid, xlab="Normal quantiles", ylab="Residuals", ...)

  # Residuals vs fitted values
  trellis.obj[["ResidFitted"]] <- xyplot(Resid ~ Fitted,
      panel = function(x, y, ...) {
        panel.xyplot(x, y, ..., type = c("p", "g"))
        panel.abline(h = 0)
        panel.loess(x, y, col="red", ...)
      }, ylab = "Residuals", xlab = "Fitted values", ...)

  # Residuals vs fitted values
  trellis.obj[["AbsResidFitted"]] <- xyplot(abs(Resid) ~ Fitted,
      panel = function(x, y, ...) {
        panel.xyplot(x, y, ..., type = c("p", "g"))
        panel.loess(x, y, col="red", ...)
      }, ylab = "|Residuals|", xlab = "Fitted values", ...)

  adt <- trellis.par.get("add.text")
  xlb <- trellis.par.get("par.xlab.text")
  ylb <- trellis.par.get("par.ylab.text")
  zlb <- trellis.par.get("par.zlab.text")
  axt <- trellis.par.get("axis.text")
  syx <- trellis.par.get("plot.symbol")
  trellis.par.set("add.text", list(cex = 0.75))
  trellis.par.set("par.xlab.text", list(cex = 0.75))
  trellis.par.set("par.ylab.text", list(cex = 0.75))
  trellis.par.set("par.zlab.text", list(cex = 0.75))
  trellis.par.set("axis.text", list(cex = 0.75))
  trellis.par.set("plot.symbol", list(cex = 0.6))

  print(trellis.obj[["histogram"]], position = c(0, 0.5, 0.5, 1), more = TRUE)
  print(trellis.obj[["qq"]], position = c(0.5, 0.5, 1, 1), more = TRUE)
  suppressWarnings(print(trellis.obj[["ResidFitted"]], position = c(0, 0, 0.5, 0.5), more = TRUE))
  suppressWarnings(print(trellis.obj[["AbsResidFitted"]], position = c(0.5, 0, 1, 0.5)))

  trellis.par.set("add.text", adt)
  trellis.par.set("par.xlab.text", xlb)
  trellis.par.set("par.ylab.text", ylb)
  trellis.par.set("par.zlab.text", zlb)
  trellis.par.set("axis.text", axt)
  trellis.par.set("plot.symbol", syx)
  invisible(trellis.obj)
}
