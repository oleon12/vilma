# Unified interactive viewer for *Vilma* spatial phylogenetic outputs

`view.vilma()` is a convenience wrapper that automatically detects the
type of *Vilma* analysis result and opens the corresponding interactive
visualization. This function removes the need to remember individual
viewer functions for each result type, providing a streamlined and
intuitive way to explore biodiversity patterns across space.

It dispatches to the appropriate map viewer based on object class:

|  |  |  |
|----|----|----|
| Class | Routed Viewer Function | Visualizes |
| `vilma.dist` | [`view.vilma.dist()`](https://oleon12.github.io/vilma/reference/view.vilma.dist.md) | Species richness & abundance rasters |
| `vilma.pd` | [`view.vilma.pd()`](https://oleon12.github.io/vilma/reference/view.vilma.pd.md) | Alpha–diversity rasters (PD, MPD, MNTD, PE, RaoQ, etc.) |
| `vilma.beta` | [`view.vilma.beta()`](https://oleon12.github.io/vilma/reference/view.vilma.beta.md) | Beta–diversity rasters (UniFrac, PhyloSor, βMPD, βMNTD, Rao β, etc.) |
| `vilma.null` | [`view.vilma.null()`](https://oleon12.github.io/vilma/reference/view.vilma.null.md) | Null–model SES maps and significance patterns |

All viewers use **leaflet** to render maps with:

- Viridis color palettes

- Pixel value pop-ups on cursor hover

- Toggleable layers

- Background basemap controls

## Usage

``` r
view.vilma(vilma)
```

## Arguments

- vilma:

  A `vilma.*` object — one of: `vilma.dist`, `vilma.pd`, `vilma.beta`,
  or `vilma.null`.

## Value

A `leaflet` HTML widget containing the appropriate interactive map. If
the object class is not recognized, an informative error is returned.

## Details

This function is especially useful during exploratory phases, teaching,
and demonstration of phylogenetic spatial patterns, allowing users to
seamlessly jump between object types without manually specifying
viewers.

## See also

`view.vilma.dist`, `view.vilma.pd`, `view.vilma.beta`, `view.vilma.null`

## Author

Omar Daniel Leon-Alvarado (<https://leon-alvarado.weebly.com>) J. Angel
Soto-Centeno (<https://www.mormoops.com>)

## Examples

``` r
if (FALSE) { # \dontrun{
# dist <- points_to_raster(occ_data)
# pd   <- faith.pd(tree, dist)
# beta <- unifrac(tree, dist)
# null <- mpd.calc.null(pd, tree, dist, iterations = 999)

view.vilma(dist)
view.vilma(pd)
view.vilma(beta)
view.vilma(null)
} # }
```
