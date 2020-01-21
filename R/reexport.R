#' @importFrom statgenSTA createTD
#' @export
statgenSTA::createTD

#' @importFrom statgenSTA fitTD
#' @export
statgenSTA::fitTD

#' @importFrom statgenSTA STAtoTD
#' @export
statgenSTA::STAtoTD

#' @importFrom utils getFromNamespace
#' @keywords internal
createSTA <- getFromNamespace(x = "createSTA", ns = "statgenSTA")

#' @keywords internal
createReport <- getFromNamespace(x = "createReport", ns = "statgenSTA")
