#' View \code{vilma.region} Rasters in an Interactive Leaflet Map
#'
#' Renders all rasters stored in a \code{vilma.region} object as toggleable
#' overlays in an interactive Leaflet map. Color scales are generated with
#' \pkg{viridisLite} and per-layer legends are added automatically. Raster
#' values can be queried on hover via \pkg{leafem}.
#'
#' @param region A \code{vilma.region} object containing a named list
#'   \code{region$raster} of \pkg{terra} \code{SpatRaster} layers (e.g.,
#'   \code{group.raster} and \code{mean.group.raster}). The function uses
#'   \code{region$raster$group.raster} to compute the initial map center when
#'   available; otherwise, it falls back to the first valid raster layer.
#'
#' @details
#' \strong{Verification:} The function checks that \code{region} is of class
#' \code{"vilma.region"} and that \code{region$raster} exists and is non-empty.
#'
#' \strong{Centering:} The initial map view is centered on the midpoint of the
#' extent of \code{region$raster$group.raster} when available. If that layer is
#' missing or invalid, the first raster layer with a valid extent is used.
#'
#' \strong{Color mapping:} For each raster, a numeric color function is created
#' with \code{viridisLite::viridis(256, option = "D")} mapped over the finite
#' value range (NAs ignored). Legends are added with the raster name as the title.
#'
#' \strong{Interaction:} Raster values are queried on mouse move using
#' \code{leafem::addImageQuery(project = TRUE)}, displaying a live readout
#' (rounded to 2 decimals) in the top-right corner.
#'
#' \strong{Layers:} Two base maps labeled "Esri" and "Carto" are added, and a
#' layers control allows toggling raster overlays. If multiple overlays are
#' available, the function hides all but the first overlay at load.
#'
#' @return A \pkg{leaflet} \code{htmlwidget} map.
#'
#' @section Performance:
#' Very large rasters can be slow in browsers. Consider aggregating
#' (\code{terra::aggregate()}) or reducing resolution prior to viewing.
#'
#' @examples
#' \dontrun{
#' reg <- vilma.regionalize(beta)
#' view.vilma.region(reg)
#' }
#'
#' @seealso
#' \code{\link[leaflet]{leaflet}},
#' \code{\link[leaflet]{addRasterImage}},
#' \code{\link[leafem]{addImageQuery}},
#' \code{\link[viridisLite]{viridis}},
#' \code{\link[terra]{ext}},
#' \code{\link{vilma.regionalize}}
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export

view.vilma.region <- function(region){

  if(!inherits(region, "vilma.region")){
    stop("region is not a vilma.region object")
  }
  
  if(is.null(region$raster) || length(region$raster) == 0){
    stop("No raster layers found")
  }
  
  ################################################################
  #                     Preparing Rasters                         #
  ################################################################
  
  ############## Helpers ###############################
  .finvals <- function(r) {
    v <- try(terra::values(r), silent = TRUE)
    if (inherits(v, "try-error") || is.null(v)) return(numeric(0))
    v[is.finite(v)]
  }
  
  # Choose center: prefer group.raster if available and valid, else first valid layer
  r_center <- NULL
  if (!is.null(region$raster$group.raster)) {
    ex_try <- try(terra::ext(region$raster$group.raster), silent = TRUE)
    if (inherits(ex_try, "SpatExtent")) r_center <- region$raster$group.raster
  }
  if (is.null(r_center)) {
    for (r in region$raster) {
      ex_try <- try(terra::ext(r), silent = TRUE)
      if (inherits(ex_try, "SpatExtent")) { r_center <- r; break }
    }
  }
  if (is.null(r_center)) stop("No valid raster extent found in region$raster")
  
  ############## Find raster center ######################
  ex <- terra::ext(r_center)
  center.x <- as.vector((ex$xmin + ex$xmax) / 2)
  center.y <- as.vector((ex$ymin + ex$ymax) / 2) 
  
  ############### Set palettes by Raster #################
  cols <- vector("list", length(region$raster))
  vals_list <- vector("list", length(region$raster))
  for (i in seq_along(region$raster)) {
    r_prj <- project(region$raster[[i]], "EPSG:4326")
    vals <- .finvals(r_prj)
    vals_list[[i]] <- vals
    
    if(length(vals)){
      rng <- range(vals, na.rm = TRUE)
      eps <- (rng[2] - rng[1]) * 1e-8
      if (!is.finite(eps) || eps == 0) eps <- 1e-8
  
      cols[[i]] <- leaflet::colorNumeric(palette  = viridisLite::viridis(256, option = "D"),
                                         domain   = c(rng[1] - eps, rng[2] + eps),
                                         na.color = "transparent",
                                         reverse  = FALSE)
    } else {
      cols[[i]] <- NULL
    }
  }
  
  # Ensure names
  nm <- names(region$raster)
  if (is.null(nm) || any(!nzchar(nm))) {
    nm <- paste0("layer_", seq_along(region$raster))
    names(region$raster) <- nm
  }
  names(cols) <- names(region$raster)
  
  ################# Base Map ############################
  map <- leaflet::leaflet(options = leaflet::leafletOptions(minZoom = 3)) %>%
    leaflet::addProviderTiles("Esri",    group = "Esri")  %>%
    leaflet::addProviderTiles("CartoDB", group = "Carto") %>%
    leaflet::setView(lng = center.x, lat = center.y, zoom = 3)
    
   ################## Add Rasters ########################
  shown_groups <- character(0)
  for (i in seq_along(cols)) {
    rname <- names(cols)[i]
    pal   <- cols[[i]]
    if (is.null(pal) || length(vals_list[[i]]) == 0) {
      map <- map %>%
        leaflet::addControl(
          html = sprintf(
            "<div style='background:#fff;padding:6px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'><b>%s</b>: no finite values (all NA)</div>",
            htmltools::htmlEscape(rname)
          ),
          position = "topright"
        )
      next
    }
    
    map <- map %>%
      leaflet::addRasterImage(x = region$raster[[i]], group = rname, colors = pal) %>%
      leaflet::addLegend(pal = pal, values = vals_list[[i]],
                         position = "bottomleft", group = rname,
                         opacity = 1, title = rname) %>%
      leafem::addImageQuery(x = region$raster[[i]], layerId = rname,
                            type = "mousemove", project = TRUE, digits = 2,
                            position = "topright", prefix = "Value: ")
    shown_groups <- c(shown_groups, rname)
  }
  
  ################# Control Layers ######################
  if (length(shown_groups)) {
    map <- map %>%
      leaflet::addLayersControl(
        baseGroups    = c("Esri","Carto"),
        overlayGroups = shown_groups
      )
    if (length(shown_groups) > 1) {
      map <- map %>% leaflet::hideGroup(setdiff(shown_groups, shown_groups[1]))
    }
  } else {
    map <- map %>%
      leaflet::addControl(
        html = "<div style='background:#fff;padding:8px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.2)'>No layers with finite values to display.</div>",
        position = "topright"
      )
  }
  
  map  


}
