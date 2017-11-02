#' Draw a biplot on the current graphics device
#'
#' Draws a biplot on the current graphics device.
#' @param x A data frame or a matrix containing the data.
#' @param centre  A logical value indicating whether the variables
#'  should be shifted to be zero centred. Alternately, a vector of
#'  length equal to the number of columns of \code{x} can be supplied.
#'  The value is passed to \code{scale}.
#' @param scale A logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes
#'  place. The default is \code{FALSE} for consistency with \code{\link{biplot.default}}, but
#'  in general scaling is advisable.  Alternatively, a vector of length
#'  equal to the number of columns of \code{x} can be supplied.  The
#'  value is passed to \code{\link{scale}}.
#' @param ... Additional arguments to be passed to \code{\link{biplot.default}}.
#' 
#' @export

ST.biplot=function(x, centre=TRUE, scale=FALSE,...){
    
    X=as.matrix(x)
    gnames=rownames(X)
    if (is.null(gnames)) gnames=1:nrow(X)
    ngnames=sapply(as.character(gnames),nchar)
    if (mean(ngnames)<=5){
      rownames(X) = gnames
    }else{
      rownames(X) = 1:length(gnames)
    }
    pca=prcomp(x=na.omit(X),center=centre, scale.=scale)
    cump2=summary(pca)$importance[3,2]
    prop.pc1=summary(pca)$importance[2,1]
    prop.pc2=summary(pca)$importance[2,2]
  
    biplot(pca, main=paste("Principal components biplot (",round(cump2*100,1),"%)",sep=""),
      xlab=paste("PC1 (",round(prop.pc1*100,1),"%)",sep=""),
      ylab=paste("PC2 (",round(prop.pc2*100,1),"%)",sep=""),...)
    
}
