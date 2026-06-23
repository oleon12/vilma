# Print summary for a `vilma.beta` object

Displays a concise summary of the results from
[`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md)
objects, including general dataset information, summary statistics for
the PhyloSor similarity and dissimilarity matrices, and ordination
diagnostics (PCoA and NMDS).

## Usage

``` r
# S3 method for class 'vilma.beta'
print(x, ...)
```

## Arguments

- x:

  An object of class `vilma.beta`, typically produced by
  [`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md).

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the input object `x`, allowing further use in
pipelines.

## Details

The print method provides an overview of the main components contained
in a `vilma.beta` object:

- The number of analyzed communities (cells) and pairwise comparisons.

- Summary statistics (mean +/- SD) for PhyloSor similarity and
  dissimilarity.

- The proportion of variance explained by the first two PCoA axes (based
  on the dissimilarity matrix).

- The NMDS stress value, indicating ordination fit quality.

- A list of available raster outputs (mean similarity/dissimilarity per
  cell, PCoA axes 1-2, and NMDS axes 1-2).

This function is automatically invoked when printing a `vilma.beta`
object to the console.

## See also

[`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md),
[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`pd.taxon`](https://oleon12.github.io/vilma/reference/pd.taxon.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)

phs <- phylosor.calc(tree, dist)
print(phs)

# The summary will display overall similarity, dissimilarity,
# and ordination statistics.
} # }
```
