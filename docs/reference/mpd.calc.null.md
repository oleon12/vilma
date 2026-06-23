# Null model for Mean Pairwise Distance (MPD) phylogenetic diversity

This function implements null models for Mean Pairwise Distance (MPD)
using different randomization approaches. It allows users to test
whether observed MPD values are higher or lower than expected under a
selected null model, either at the global level across all cells or at
the cell level.

## Usage

``` r
mpd.calc.null(
  pd,
  tree,
  dist,
  iterations = 999,
  method = c("global", "cell"),
  sampling = c("taxa.label", "range", "neighbor", "regional"),
  n.directions = c("rook", "bishop", "queen"),
  regional.weight = c("uniform", "frequency", "range")
)
```

## Arguments

- pd:

  An object of class `vilma.pd`, containing observed MPD values
  calculated with
  [`mpd.calc`](https://oleon12.github.io/vilma/reference/mpd.calc.md).

- tree:

  A rooted phylogenetic tree of class `phylo`, with branch lengths.

- dist:

  An object of class `vilma.dist`, representing species distributions
  across spatial cells. See
  [`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md),
  [`rast2occ`](https://oleon12.github.io/vilma/reference/rast2occ.md),
  or [`pol2occ`](https://oleon12.github.io/vilma/reference/pol2occ.md).

- iterations:

  Integer. Number of randomizations to perform. The default is `999`.

- method:

  Character. Null-model output method. Options are:

  - `"global"`: evaluates the total observed MPD across all cells
    against a null distribution of total MPD values.

  - `"cell"`: evaluates observed MPD independently for each cell and
    returns cell-level standardized effect sizes and p-values.

- sampling:

  Character. Randomization strategy. Options are:

  - `"taxa.label"`: randomly permutes species labels on the phylogenetic
    tree while keeping the spatial distribution fixed.

  - `"range"`: randomizes the species-by-cell matrix using a swap
    algorithm while preserving matrix structure.

  - `"neighbor"`: swaps species occurrences between adjacent spatial
    cells.

  - `"regional"`: samples species from the regional pool while
    preserving observed cell richness.

- n.directions:

  Character. Neighborhood adjacency definition used only when
  `sampling = "neighbor"`. Options are `"rook"`, `"bishop"`, and
  `"queen"`. The default is `"queen"`.

- regional.weight:

  Character. Weighting scheme for the regional species pool used only
  when `sampling = "regional"`. Options are:

  - `"uniform"`: all species have equal sampling probability.

  - `"frequency"`: species are weighted by their number of occurrence
    records.

  - `"range"`: species are weighted by the number of occupied cells.

  The default is `"uniform"`.

## Value

An object of class `vilma.null`. The returned object depends on the
selected `method`.

For `method = "global"`, the object contains:

- `pd.obs`: observed total MPD across all cells.

- `null.pd`: null distribution of total MPD values.

- `SES`: standardized effect size.

- `Pvalue`: p-value calculated from the null distribution.

- `Iterations`: number of randomizations.

- `Iter.table`: table of MPD values for each iteration.

- `Method`: selected null-model method.

For `method = "cell"`, the object contains:

- `CellValues`: data frame with observed MPD, null mean, null standard
  deviation, SES, and p-value for each cell.

- `Raster`: raster layer containing cell-level SES values.

- `Iterations`: number of randomizations.

- `Iter.table`: table of MPD values for each iteration.

- `Method`: selected null-model method.

## Details

The function generates a null distribution of MPD values based on the
selected randomization strategy.

When `method = "global"`, the observed MPD values are summed across all
cells and compared with a null distribution generated from the sum of
MPD values across randomized iterations. The function returns the
observed total MPD, the null distribution, a standardized effect size,
and a p-value.

When `method = "cell"`, the function compares the observed MPD of each
cell against a cell-specific null distribution. It returns observed MPD,
null mean, null standard deviation, standardized effect size, p-value,
and a raster containing the cell-level SES values.

For `sampling = "taxa.label"`, species labels are randomized on the
phylogeny, while the spatial distribution remains fixed. This approach
evaluates whether the observed association between species identities
and phylogenetic distances differs from random expectations.

For `sampling = "range"`, the species-by-cell matrix is randomized using
a swap procedure. This approach changes species occurrences among cells
while preserving key properties of the incidence matrix.

For `sampling = "neighbor"`, species are swapped among adjacent cells.
The neighborhood structure is defined by `n.directions`, allowing rook,
bishop, or queen adjacency. This approach is useful when spatially
constrained randomizations are desired.

For `sampling = "regional"`, species are sampled from the regional
species pool while preserving the observed richness of each cell.
Species can be sampled with equal probability or weighted by frequency
or range size.

## Interpretation

Negative SES values indicate phylogenetic clustering, meaning that
observed MPD is lower than expected under the selected null model.
Positive SES values indicate phylogenetic overdispersion, meaning that
observed MPD is higher than expected. SES values near zero suggest that
observed MPD is close to the null expectation.

The interpretation of SES and p-values depends on the selected
randomization strategy. For example, a taxa-label null model tests
whether the observed phylogenetic structure differs from random
associations between species names and the tree, whereas range,
neighbor, and regional-pool models evaluate different assumptions about
spatial occurrence structure.

## References

Gotelli, N. J., & Graves, G. R. (1996). *Null Models in Ecology*.
Smithsonian Institution Press.

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. *Ecology*, 81(9), 2606–2621.
[doi:10.1890/0012-9658(2000)081\[2606:NMAOSC\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282000%29081%5B2606%3ANMAOSC%5D2.0.CO%3B2)

Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002).
Phylogenies and community ecology. *Annual Review of Ecology and
Systematics*, 33, 475–505.
[doi:10.1146/annurev.ecolsys.33.010802.150448](https://doi.org/10.1146/annurev.ecolsys.33.010802.150448)

## See also

[`mpd.calc`](https://oleon12.github.io/vilma/reference/mpd.calc.md),
[`mntd.calc.null`](https://oleon12.github.io/vilma/reference/mntd.calc.null.md),
[`faith.pd.null`](https://oleon12.github.io/vilma/reference/faith.pd.null.md),
[`pe.calc.null`](https://oleon12.github.io/vilma/reference/pe.calc.null.md),
[`rao.calc.null`](https://oleon12.github.io/vilma/reference/rao.calc.null.md),
[`swap.null`](https://oleon12.github.io/vilma/reference/swap.null.md),
[`return.vilma.dist`](https://oleon12.github.io/vilma/reference/return.vilma.dist.md),
[`return.vilma.dist2`](https://oleon12.github.io/vilma/reference/return.vilma.dist2.md),
[`view.vilma.null`](https://oleon12.github.io/vilma/reference/view.vilma.null.md),
[`plot.vilma.null`](https://oleon12.github.io/vilma/reference/plot.vilma.null.md),
[`write.vilma.null`](https://oleon12.github.io/vilma/reference/write.vilma.null.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/>

J. Angel Soto-Centeno <https://www.mormoops.com/>

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_tree)
data(example_dist)

pd <- mpd.calc(
  tree = example_tree,
  dist = example_dist
)

null_model <- mpd.calc.null(
  pd = pd,
  tree = example_tree,
  dist = example_dist,
  iterations = 99,
  method = "global",
  sampling = "taxa.label"
)

print(null_model)
plot(null_model)
view.vilma(null_model)
} # }
```
