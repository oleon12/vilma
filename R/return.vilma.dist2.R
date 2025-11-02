#' Convert a VILMA Distribution Matrix to a Data Frame
#' @description
#' Converts a list of presence per cell into a long-format
#' data frame compatible with VILMA functions. Each row represents a species-cell record.
#'
#' @param x A lisy object representing species distributions (presence).
#'
#' @details
#' This function takes a matrix where rows are cells and columns are species,
#' and returns a data frame with columns \code{Sp}, \code{Lon}, \code{Lat}, and \code{Cell}.
#' The \code{Lon} and \code{Lat} columns are filled with zeros by default, 
#' as coordinates are not required for internal VILMA calculations.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{Sp}{Species name or identifier.}
#'   \item{Lon}{Longitude (default 0).}
#'   \item{Lat}{Latitude (default 0).}
#'   \item{Cell}{Cell identifier corresponding to the row in the input matrix.}
#' }
#'
#' @author 
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export

return.vilma.dist2 <- function(x){
  
  out <- data.frame(Sp = NA,
                    Lon = NA,
                    Lat = NA,
                    Cell = NA)
  
  for(i in 1:length(x)){
    y <- data.frame(Sp = x[[i]],
                    Lon = rep(0, length(x[[i]])),
                    Lat = rep(0, length(x[[i]])),
                    Cell = rep(names(x)[i], length(x[[i]])))
   
   out <- rbind(out, y)
   
  }
  
  out <- out[-1,]
  
  return(out)


}
