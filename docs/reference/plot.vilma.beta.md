# Plot method for `vilma.beta` objects

Produces quick-look maps of mean PhyloSor similarity and mean PhyloSor
dissimilarity per cell from a `vilma.beta` object (output of
[`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md)).
Values are drawn using the underlying `terra` raster and annotated with
the cell-wise means.

## Usage

``` r
# S3 method for class 'vilma.beta'
plot(x, ...)
```

## Arguments

- x:

  An object of class `vilma.beta`, typically returned by
  [`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md)
  and containing the rasters `$rasters$mean.similarity` and
  `$rasters$mean.dissimilarity`.

- ...:

  Additional graphical arguments (currently ignored).

## Value

Invisibly returns `x` (the input object), allowing further use in
pipelines.

## Details

The function opens a 2-panel plotting layout (`par(mfrow=c(2,1))`) and
renders:

1.  Mean PhyloSor similarity per cell (`$rasters$mean.similarity`).

2.  Mean PhyloSor dissimilarity per cell
    (`$rasters$mean.dissimilarity`).

For each panel, numeric values are overlaid at cell centroids (rounded
to 2 decimals). The plotting is performed with terra (i.e.,
[`terra::plot`](https://rspatial.github.io/terra/reference/plot.html)).

Expect `NA` cells when communities were filtered (e.g.,
`method = "exclude"`) or when no value could be computed. The function
temporarily modifies graphical parameters and restores them on exit.

## See also

[`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md),
[`print.vilma.beta`](https://oleon12.github.io/vilma/reference/print.vilma.beta.md),
[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>
