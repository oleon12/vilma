#' Convert a VILMA Distribution Matrix to a Data Frame
#' @description
#' Converts a presence-absence matrix (or similar object) into a long-format
#' data frame compatible with VILMA functions. Each row represents a species-cell record.
#'
#' @param x A matrix or object representing species distributions (presence/absence).
#' @param ... Additional arguments (currently unused).
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
#' @author Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export

return.vilma.dist <- function(x, ...){
  #Cells0 <- rownames(x)
  
  #Cells <- c()
  #Sp <- c()
  
  #for(i in 1:length(Cells0)){
  #  Cells <- c(Cells, rep(Cells0[i], length(which(x[i,]==1))))
  #  Sp <- c(Sp, colnames(x)[which(x[i,]==1)])
  #}
  
  #out <- data.frame(Sp= Sp,
  #                  Lon = rep(0, length(Sp)),
  #                  Lat = rep(0, length(Sp)),
  #                  Cell = Cells)
  

df <- as.data.frame(as.table(x))
df <- df[df$Freq > 0, ]
df <- data.frame(Sp = df$Var2,
                 Lon = rep(0, nrow(df)),
                 Lat = rep(0, nrow(df)),
                 Cell = df$Var1)
return(df)
  
  
  #return(out)
  
}
