# Net Relatedness Index (NRI)

Computes the Net Relatedness Index (NRI) for communities (grid cells)
based on phylogenetic structure. NRI is derived from the standardized
effect size (SES) of the Mean Pairwise Distance (MPD) under a null
model, with the sign inverted so that positive values indicate
phylogenetic clustering and negative values indicate phylogenetic
overdispersion.

## Usage

``` r
nri.calc(
  tree,
  dist,
  mpd.method = c("root", "node", "exclude"),
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

- mpd.method:

  Character string specifying how to handle single-species cells in MPD.
  Options are:

  - `"root"` – For cells with single species, uses distance from root to
    taxon.

  - `"node"` – For cells with single species, uses branch length to
    nearest ancestral node.

  - `"exclude"` – Excludes single-species cells from calculations
    (default).

- abundance:

  Logical. If `TRUE`, weights MPD by relative species abundances within
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

- `pd.table` – Data frame with cells, species richness, and NRI values.

- `rasters` – List with abundance, richness, and NRI rasters.

- `calculation.method` – Method used for single-species cells in MPD.

- `abundance` – Logical, whether abundance-weighting was applied.

- `index` – The string `"nri.calc"`.

## Details

The Net Relatedness Index (NRI) is defined as: \$\$NRI = -
\frac{MPD\_{obs} - \mu(MPD\_{null})}{\sigma(MPD\_{null})}\$\$ Positive
NRI values indicate that species within communities are more closely
related (phylogenetic clustering) than expected under the null model.
Negative NRI values indicate that species are less related (phylogenetic
overdispersion).

NRI is computed at the `"cell"` level, returning results per grid cell,
including a raster of NRI values for visualization.

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

[`mpd.calc`](https://oleon12.github.io/vilma/reference/mpd.calc.md),
[`mpd.calc.null`](https://oleon12.github.io/vilma/reference/mpd.calc.null.md),
[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com/>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)

# Calculate NRI using default null model (taxa.label)
nri_result <- nri.calc(tree = tree, dist = dist, iterations = 999)

print(nri_result)
plot(nri_result)
view.vilma(nri_result)

# Calculate NRI using neighbor null model
nri_neighbor <- nri.calc(tree = tree, dist = dist,
                         iterations = 999, sampling = "neighbor",
                         n.directions = "queen")
} # }


```
