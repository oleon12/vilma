#' Plot VILMA Correlation Results
#'
#' Provides a quick visualization of a \code{vilma.cor} object by plotting the
#' relationship between species richness and phylogenetic diversity (PD) values.
#' The method automatically selects the best-supported model based on the minimum
#' AIC value stored in \code{x$AIC}. If \pkg{ggplot2} is available, a scatterplot
#' with a fitted trend line is produced; otherwise, the function falls back to a
#' base R plot using the selected fitted model object.
#'
#' @param x An object of class \code{vilma.cor}.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' Model selection is performed by choosing the model name corresponding to the
#' smallest value in \code{x$AIC}. When \pkg{ggplot2} is installed, the function
#' plots \code{x$Data} using \code{sr.vals} on the x-axis and \code{pd.vals} on the
#' y-axis, and overlays a fitted smooth with \code{geom_smooth(method = <best model>)}.
#' If \pkg{ggplot2} is not installed, the function plots the best model stored in
#' \code{x$Models[[method]]} using base R.
#'
#' @return A plot is produced. Invisibly returns \code{x}.
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com/}
#'
#' @seealso \code{\link{vilma.cor}}
#'
#' @export
#' @method plot vilma.cor

plot.vilma.cor <- function(x, ...) {

  if(!inherits(x, "vilma.cor")){
    stop("Input must be a 'vilma.cor' object")
  }

  method <- as.character(names(x$AIC)[which(x$AIC == min(x$AIC))])

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(x$Data) +
      ggplot2::geom_point(ggplot2::aes(x = .data$sr.vals, y = .data$pd.vals)) +
      ggplot2::xlab("Species Richness") +
      ggplot2::ylab("PD values") +
      ggplot2::geom_smooth(ggplot2::aes(x = .data$sr.vals, y = .data$pd.vals),
                           method = method) +
      ggplot2::theme_bw()
  } else {
    graphics::plot(x$Models[[method]], xlab="Species Richness", ylab="PD values")
  }

  invisible(x)
}
