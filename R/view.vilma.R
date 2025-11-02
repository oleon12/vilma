#' Unified interactive viewer for *Vilma* spatial phylogenetic outputs
#'
#' @description
#' `view.vilma()` is a convenience wrapper that automatically detects the type of
#' *Vilma* analysis result and opens the corresponding interactive visualization.
#' This function removes the need to remember individual viewer functions for
#' each result type, providing a streamlined and intuitive way to explore
#' biodiversity patterns across space.
#'
#' It dispatches to the appropriate map viewer based on object class:
#'
#' | Class          | Routed Viewer Function         | Visualizes |
#' |----------------|--------------------------------|-----------|
#' | `vilma.dist`   | `view.vilma.dist()`            | Species richness & abundance rasters |
#' | `vilma.pd`     | `view.vilma.pd()`              | Alpha–diversity rasters (PD, MPD, MNTD, PE, RaoQ, etc.) |
#' | `vilma.beta`   | `view.vilma.beta()`            | Beta–diversity rasters (UniFrac, PhyloSor, βMPD, βMNTD, Rao β, etc.) |
#' | `vilma.null`   | `view.vilma.null()`            | Null–model SES maps and significance patterns |
#'
#' All viewers use **leaflet** to render maps with:
#' - Viridis color palettes
#' - Pixel value pop-ups on cursor hover
#' - Toggleable layers
#' - Background basemap controls
#'
#' @param vilma
#' A `vilma.*` object — one of:  
#' `vilma.dist`, `vilma.pd`, `vilma.beta`, or `vilma.null`.
#'
#' @return
#' A `leaflet` HTML widget containing the appropriate interactive map.
#' If the object class is not recognized, an informative error is returned.
#'
#' @details
#' This function is especially useful during exploratory phases, teaching,
#' and demonstration of phylogenetic spatial patterns, allowing users to
#' seamlessly jump between object types without manually specifying viewers.
#'
#' @examples
#' \dontrun{
#' # dist <- points_to_raster(occ_data)
#' # pd   <- faith.pd(tree, dist)
#' # beta <- unifrac(tree, dist)
#' # null <- mpd.calc.null(pd, tree, dist, iterations = 999)
#'
#' view.vilma(dist)
#' view.vilma(pd)
#' view.vilma(beta)
#' view.vilma(null)
#' }
#'
#' @seealso
#' `view.vilma.dist`, `view.vilma.pd`, `view.vilma.beta`, `view.vilma.null`
#'
#' @author
#' Omar Daniel Leon-Alvarado (<https://leon-alvarado.weebly.com>)  
#' J. Angel Soto-Centeno (<https://www.mormoops.com>)
#'
#' @family Vilma visualization functions
#'
#' @export

view.vilma <- function(vilma){

  if(inherits(vilma, "vilma.dist")){
  
    map <- view.vilma.dist(vilma)
  
  }
  
  if(inherits(vilma, "vilma.pd")){
  
    map <- view.vilma.pd(vilma)
  
  }
  
  if(inherits(vilma, "vilma.beta")){
  
    map <- view.vilma.beta(vilma)
  
  }
  
  if(inherits(vilma, "vilma.null")){
  
    map <- view.vilma.null(vilma)
  
  }
  
  map

}
