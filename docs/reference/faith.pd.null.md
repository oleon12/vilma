# Null model for Faith's Phylogenetic Diversity (PD)

This function implements null models for Faith's PD using different
randomization approaches. It allows testing significance of PD at either
the global (across all cells) or cell level, with several sampling
strategies: randomizing taxa labels, species ranges, or neighborhood
swaps.

## Usage

``` r
faith.pd.null(
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

  An object of class `vilma.pd`, containing observed PD values.

- tree:

  A rooted phylogenetic tree of class `phylo`, with branch lengths.

- dist:

  An object of class `vilma.dist`, representing species distributions
  (see `points.to.raster()`).

- iterations:

  Integer. Number of randomizations to perform (default = 999).

- method:

  Character. Null model output method:

  - `"global"` – significance of total PD across all cells.

  - `"cell"` – significance per cell.

- sampling:

  Character. Randomization strategy:

  - `"taxa.label"` – permutes species labels on the tree.

  - `"range"` – randomizes species ranges using swap algorithms.

  - `"neighbor"` – swaps species occurrences between adjacent cells.

- n.directions:

  Character. Neighborhood adjacency definition for the `"neighbor"`
  method. Options: `"rook"`, `"bishop"`, `"queen"` (default =
  `"queen"`).

- regional.weight:

  Weighting of the regional pool for the `"regional"` model: `"uniform"`
  (equal), `"frequency"` (by number of records), or `"range"` (by number
  of occupied cells). Defaults to `"uniform"`.

## Value

An object of class `vilma.null`, with structure depending on `method`:

- For `global`: a list with observed PD, null distribution, SES,
  p-value, iteration table, and metadata.

- For `cell`: a list with per-cell values (PD, null mean, null SD, SES,
  p-value), a raster of SES values, iteration table, and metadata.

## Details

The function generates a null distribution of PD values based on the
chosen sampling strategy. For the `"global"` method, the observed PD is
compared against the distribution of total PD values across iterations.
For the `"cell"` method, standardized effect sizes (SES) and p-values
are computed for each grid cell individually.

## References

Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity.
*Biological Conservation*, 61(1), 1–10.
[doi:10.1016/0006-3207(92)91201-3](https://doi.org/10.1016/0006-3207%2892%2991201-3)

Gotelli, N. J., & Graves, G. R. (1996). *Null Models in Ecology*.
Smithsonian Institution Press.

Gotelli, N. J. (2000). Null model analysis of species co‐occurrence
patterns. *Ecology*, 81(9), 2606–2621.
[doi:10.1890/0012-9658(2000)081\[2606:NMAOSC\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282000%29081%5B2606%3ANMAOSC%5D2.0.CO%3B2)

Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002).
Phylogenies and community ecology. *Annual Review of Ecology and
Systematics*, 33, 475–505.
[doi:10.1146/annurev.ecolsys.33.010802.150448](https://doi.org/10.1146/annurev.ecolsys.33.010802.150448)

## See also

[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md),

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com/>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)
pd <- faith.pd(tree, dist, method = "cell")
null_model <- faith.pd.null(pd, tree, dist,
                            iterations = 999, method = "cell",
                            sampling = "range")
} # }

```
