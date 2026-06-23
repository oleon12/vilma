# Interactive viewer for vilma null-model results (SES raster)

Renders a **leaflet** map to explore the SpatialRaster of standardized
effect sizes (SES) produced by a `vilma.null` object (e.g., from
[`rao.calc.null`](https://oleon12.github.io/vilma/reference/rao.calc.null.md)).
The map centers on the raster extent, uses a `viridis` palette, shows a
live value readout on mouse move, and provides base-layer and legend
controls.

## Usage

``` r
view.vilma.null(null)
```

## Arguments

- null:

  A `vilma.null` object containing a SES raster in `null$Raster`. The
  function will error if `null$Raster` is missing or `NULL`.

## Value

A `leaflet` `htmlwidget` map for interactive viewing.

## Details

- Centers the map at the midpoint of `terra::ext(null$Raster)`.

- Uses `viridisLite::viridis(256, option = "D")` via
  [`leaflet::colorNumeric`](https://rstudio.github.io/leaflet/reference/colorNumeric.html)
  for continuous SES coloring (`na.color = NA`).

- Adds `"Esri"` and `"CartoDB"` provider tiles as base layers.

- Displays a legend titled *"SES values"* (bottom-left).

- Displays an on-hover pixel value panel using
  [`leafem::addImageQuery`](https://r-spatial.github.io/leafem/reference/addImageQuery.html)
  (top-right).

- Map starts at `minZoom = 3`; initial zoom is set to 3.

## See also

[`rao.calc.null`](https://oleon12.github.io/vilma/reference/rao.calc.null.md),
[`rao.calc`](https://oleon12.github.io/vilma/reference/rao.calc.md),
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- example_tree()
dist <- example_dist()
dist <- points_to_raster(dist)
out.c <- pe.calc.null(tree, dist)
view.vilma.null(out.c)
} # }
```
