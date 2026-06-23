# Interactive viewer for vilma β-diversity rasters

Renders a **leaflet** map to explore raster layers contained in a
`vilma.beta` object (e.g., mean Rao β dissimilarity, PCoA axes, nMDS
axes). The map centers on the extent of the first raster, uses a
`viridis` palette per layer, shows on-hover pixel values, and provides
base-layer and legend controls. All overlay rasters are added as
separate toggleable groups; by default only the first group is shown.

## Usage

``` r
view.vilma.beta(beta)
```

## Arguments

- beta:

  A `vilma.beta` object with a named list of raster layers in
  `beta$rasters` (e.g., `mean.dissimilarity`, `pcoa.1`, `pcoa.2`,
  `ndms.1`, `ndms.2`). The function will error if `beta$rasters` is
  `NULL`.

## Value

A `leaflet` `htmlwidget` with the β-diversity rasters and controls.

## Details

- Centers the map at the midpoint of `terra::ext(beta$rasters[[1]])`.

- For each raster, builds a continuous color function with
  `viridis(256, option = "D")` via
  [`leaflet::colorNumeric`](https://rstudio.github.io/leaflet/reference/colorNumeric.html)
  (`na.color = NA`).

- Adds `"Esri"` and `"CartoDB"` provider tiles as base layers.

- Adds a legend titled with the raster name (bottom-left).

- Shows an on-hover value panel using
  [`leafem::addImageQuery`](https://r-spatial.github.io/leafem/reference/addImageQuery.html)
  (top-right, `digits = 2`).

- Uses `minZoom = 3`; initial zoom is 3. All rasters are added as
  overlay groups; groups other than the first are hidden initially.

## See also

[`rao.beta`](https://oleon12.github.io/vilma/reference/rao.beta.md) for
computing Rao β dissimilarity,
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)
for building `vilma.dist`

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- example_tree()
dist <- example_dist()
dist <- points_to_raster(dist)
b <- bmpd(tree, dist)
view.vilma.beta(b)
} # }
```
