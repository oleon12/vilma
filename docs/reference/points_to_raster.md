# Convert Species Occurrence Points to Raster Distribution

Converts species occurrence points (Species, Longitude, Latitude) to a
`vilma.dist` object that includes a template grid and (optionally)
raster layers of species richness and abundance per cell.

## Usage

``` r
points_to_raster(
  points,
  crs = 4326,
  ext = NULL,
  res = 1,
  doRast = TRUE,
  symmetrical = FALSE
)
```

## Arguments

- points:

  A `matrix` or `data.frame` with exactly three columns, in the order:
  Species, Longitude, Latitude.

- crs:

  Coordinate Reference System. Either an EPSG code as integer (e.g.,
  `4326` for WGS84) or a PROJ string. Defaults to EPSG:4326.

- ext:

  A [`terra::ext`](https://rspatial.github.io/terra/reference/ext.html)
  extent defining output raster bounds. If `NULL` (default), the extent
  is computed from `points` and buffered slightly.

- res:

  Numeric resolution of the output grid in CRS units. Either a single
  value (square cells; `xres = yres`) or a length-2 numeric vector
  `c(xres, yres)`. Defaults to `1` (degrees if CRS is geographic).

- doRast:

  Logical; if `TRUE` (default) raster layers are created (richness and
  abundance). If `FALSE`, returns only the distribution table and the
  grid definition.

- symmetrical:

  Logical; if `TRUE`, forces square pixels by setting `yres <- xres` and
  snaps the extent so that its width/height are exact multiples of the
  (possibly forced) resolution. Default `FALSE`.

## Value

An object of class `vilma.dist` with:

- `distribution`: `data.frame` with original points and their `Cell`
  IDs.

- `grid`:
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  template grid.

- `r.raster`: richness raster (`#` unique species per cell) or `NULL`.

- `ab.raster`: abundance raster (`#` records per cell) or `NULL`.

## Details

**Resolution and symmetry:** `res` can be a scalar or `c(xres, yres)`.
When `symmetrical = TRUE`, the function enforces `yres = xres` and snaps
the extent to whole-cell multiples, ensuring symmetric (square) pixels.

**On geographic CRSs:** If the CRS is geographic (lon/lat), square cells
are *square degrees*, not equal-area squares. For equal-area, project to
an appropriate (equal-area) CRS before calling this function.

**Outputs:** Cell IDs are assigned via
[`terra::cellFromXY`](https://rspatial.github.io/terra/reference/xyCellFrom.html).
If `doRast = TRUE`, abundance is computed as record counts per cell and
richness as the number of unique species per cell.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Sample occurrence data
tree <- examplae_tree()
dist <- example_dist()

# Rectangular cells (0.5 x 0.25 degrees)
d2 <- points_to_raster(dist, res = c(0.5, 0.25))

# Force square cells with snapping (use xres = 0.25)
d3 <- points_to_raster(dist, res = 0.25, symmetrical = TRUE)

# Return only table + grid (no rasters)
d4 <- points_to_raster(dist, doRast = FALSE)
} # }
```
