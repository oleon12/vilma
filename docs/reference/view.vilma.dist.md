# Interactive viewer for vilma distributions (richness & abundance)

Renders a **leaflet** map to explore the richness and abundance rasters
stored in a `vilma.dist` object (e.g., from
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)).
The map centers on the raster extent, uses `viridis` palettes, shows
on-hover cell values, and provides base-layer and legend controls. By
default, the *Richness* layer is shown and *Abundance* is hidden.

## Usage

``` r
view.vilma.dist(dist)
```

## Arguments

- dist:

  A `vilma.dist` object containing `dist$r.raster` (richness) and
  `dist$ab.raster` (abundance). The function will error if either raster
  is missing or `NULL`.

## Value

A `leaflet` `htmlwidget` map with richness and abundance overlays.

## Details

- Centers the map at the midpoint of `terra::ext(dist$r.raster)`.

- Builds continuous color scales with
  `viridisLite::viridis(256, option = "D")` via
  [`leaflet::colorNumeric`](https://rstudio.github.io/leaflet/reference/colorNumeric.html)
  (`na.color = NA`).

- Adds `"Esri"` and `"CartoDB"` provider tiles as base layers.

- Adds legends titled *"richness"* and *"abundance"* (bottom-left).

- Shows an on-hover value panel for each raster using
  [`leafem::addImageQuery`](https://r-spatial.github.io/leafem/reference/addImageQuery.html)
  (top-right, `digits = 2`).

- Uses `minZoom = 3`; initial zoom is 3. *Abundance* layer is hidden
  initially.

## See also

[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)
to create `vilma.dist` objects,
[`view.vilma.beta`](https://oleon12.github.io/vilma/reference/view.vilma.beta.md),
[`view.vilma.null`](https://oleon12.github.io/vilma/reference/view.vilma.null.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'd' is a vilma.dist produced by points_to_raster(...)
view.vilma.dist(d)
} # }
```
