# Rao's Q (within-community) phylogenetic diversity

Computes Rao’s quadratic entropy (RaoQ) per spatial cell from a rooted
phylogenetic tree and a `vilma.dist` distribution. RaoQ is the expected
phylogenetic dissimilarity between two individuals randomly drawn from
the same community (cell), using patristic distances from the tree.

## Usage

``` r
rao.calc(tree, dist, abundance = FALSE, scale01 = TRUE)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths.

- dist:

  An object of class `vilma.dist` (e.g., from `points.to.raster()`),
  whose `$distribution` contains at least `Cell` and `Sp`, and whose
  `$grid` is used to rasterize outputs.

- abundance:

  Logical. If `FALSE` (default), a presence–absence Cell×Sp table is
  used (each species weighted equally within a cell). If `TRUE`, the
  Cell×Sp counts from `table(Cell, Sp)` are used to compute within-cell
  relative abundances.

- scale01:

  Logical. If `TRUE` (default), the patristic distance matrix is divided
  by its maximum off-diagonal so that \\D \in \[0,1\]\\ (diagonal set to
  0). Under this scaling, RaoQ lies in \\\[0,1\]\\ up to numerical
  tolerance.

## Value

An object of class `vilma.pd`:

- `distribution` – Original `dist$distribution`.

- `grid` – Original `dist$grid`.

- `pd.table` – Data frame with columns: `Cell` (cell id), `SR` (species
  richness), and `PD` (RaoQ).

- `rasters` – List with `ab.raster`, `r.raster`, and `pd.raster` (RaoQ
  per cell).

- `calculation.method` – Character summary of the settings.

- `index` – The string `"rao"`.

## Details

Species are intersected between `tree$tip.label` and
`dist$distribution$Sp`; only the common species are used. The patristic
matrix \\D\\ is computed with
[`ape::cophenetic.phylo()`](https://rdrr.io/pkg/ape/man/cophenetic.phylo.html)
and then **reindexed to match the column order of the community matrix**
to ensure correct behavior under taxa-label null models.

For each cell, let \\p\\ be the vector of within-cell relative
abundances (or presence–absence proportions). RaoQ is computed as \\Q =
p^\top D p\\. Empty cells (no species) are returned as `NA`.

## Interpretation

Larger RaoQ indicates greater within-community phylogenetic dispersion.
With `scale01 = TRUE`, values are comparable across trees (range approx.
`[0,1]`).

## References

Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified
approach. *Theoretical Population Biology*, 21, 24–43.  
Botta-Dukát, Z. (2005). Rao’s quadratic entropy as a measure of
functional diversity. *Community Ecology*, 6, 283–290.

## See also

[`rao.beta`](https://oleon12.github.io/vilma/reference/rao.beta.md) for
between-community dissimilarity;
[`rao.calc.null`](https://oleon12.github.io/vilma/reference/rao.calc.null.md)
for null-model SES;
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)
to build `vilma.dist`.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{

tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)
rao <- rao.calc(tree = tree, dist = dist, abundance = FALSE, scale01 = TRUE)
print(rao)
plot(rao)
view.vilma(rao)
} # }
```
