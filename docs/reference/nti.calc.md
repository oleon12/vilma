# Nearest Taxon Index (NTI)

Computes the Nearest Taxon Index (NTI) for communities (grid cells)
based on phylogenetic structure. NTI is derived from the standardized
effect size (SES) of the Mean Nearest Taxon Distance (MNTD) under a null
model, with the sign inverted so that positive values indicate
phylogenetic clustering near the tips of the tree and negative values
indicate phylogenetic overdispersion.

## Usage

``` r
nti.calc(
  tree,
  dist,
  mntd.method = c("root", "node", "exclude"),
  abundance = FALSE,
  iterations = 999,
  sampling = c("taxa.label", "range", "neigbor", "regional"),
  n.directions = c("rook", "bishop", "queen"),
  regional.weight = c("uniform", "frequency", "range")
)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo`, with branch lengths.

- dist:

  An object of class `vilma.dist`, typically produced by
  `points.to.raster()`, containing species distribution data.

- mntd.method:

  Character string specifying how to handle single-species cells in
  MNTD. Options are:

  - `"root"` – For cells with single species, uses distance from root to
    taxon.

  - `"node"` – For cells with single species, uses branch length to
    nearest ancestral node.

  - `"exclude"` – Excludes single-species cells from calculations
    (default).

- abundance:

  Logical. If `TRUE`, weights MNTD by relative species abundances within
  each cell. Default is `FALSE`.

- iterations:

  Integer. Number of randomizations for the null model (default = 999).

- sampling:

  Character. Randomization strategy for null model:

  - `"taxa.label"` – Permutes species labels on the tree.

  - `"range"` – Randomizes species ranges using swap algorithms.

  - `"neighbor"` – Swaps species occurrences between adjacent cells.

  - `"regional"` – Assembles communities by drawing from a regional
    pool.

- n.directions:

  Character. Neighborhood adjacency definition for `"neighbor"` method.
  Options: `"rook"`, `"bishop"`, or `"queen"` (default = `"queen"`).

- regional.weight:

  Character. Weighting scheme for species in the regional null model:

  - `"uniform"` – All species have equal probability.

  - `"frequency"` – Weighted by the number of records per species.

  - `"range"` – Weighted by the number of unique cells occupied.

## Value

An object of class `vilma.pd`, which is a list containing:

- `distribution` – Original species distribution.

- `grid` – The reference raster grid.

- `pd.table` – Data frame with cells, species richness, and NTI values.

- `rasters` – List with abundance, richness, and NTI rasters.

- `calculation.method` – Method used for single-species cells in MNTD.

- `abundance` – Logical, whether abundance-weighting was applied.

- `index` – The string `"nti.calc"`.

## Details

The Nearest Taxon Index (NTI) is defined as: \$\$NTI = -
\frac{MNTD\_{obs} - \mu(MNTD\_{null})}{\sigma(MNTD\_{null})}\$\$
Positive NTI values indicate that species within communities are more
closely related at the tips of the phylogeny (phylogenetic clustering)
than expected under the null model. Negative NTI values indicate
phylogenetic overdispersion at the tips.

NTI is computed at the `"cell"` level, returning results per grid cell,
including a raster of NTI values for visualization.

## References

Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002).
Phylogenies and community ecology. *Annual Review of Ecology and
Systematics*, 33, 475–505.
[doi:10.1146/annurev.ecolsys.33.010802.150448](https://doi.org/10.1146/annurev.ecolsys.33.010802.150448)

Gotelli, N. J., & Graves, G. R. (1996). *Null Models in Ecology*.
Smithsonian Institution Press.

Gotelli, N. J. (2000). Null model analysis of species co‐occurrence
patterns. *Ecology*, 81(9), 2606–2621.
[doi:10.1890/0012-9658(2000)081\[2606:NMAOSC\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282000%29081%5B2606%3ANMAOSC%5D2.0.CO%3B2)

## See also

[`mntd.calc`](https://oleon12.github.io/vilma/reference/mntd.calc.md),
[`mntd.calc.null`](https://oleon12.github.io/vilma/reference/mntd.calc.null.md),
[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)

# Calculate NTI using default null model (taxa.label)
nti_result <- nti.calc(tree = tree, dist = dist, iterations = 999)
print(nti_result)
plot(nti_result)
view.vilma(nti_result)

# Calculate NTI using regional null model with frequency weights
nti_regional <- nti.calc(tree = tree, dist = dist,
                         iterations = 999, sampling = "regional",
                         regional.weight = "frequency")
} # }
```
