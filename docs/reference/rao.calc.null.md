# Rao's Q (α) null models for vilma objects

Generates null expectations for **within-community** Rao’s quadratic
entropy (RaoQ; phylogenetic α-diversity) using several randomization
schemes, and computes standardized effect sizes (SES) either
**globally** or **per cell**.

## Usage

``` r
rao.calc.null(
  pd,
  tree,
  dist,
  iterations = 999,
  method = c("global", "cell"),
  sampling = c("taxa.label", "range", "neighbor", "regional"),
  n.directions = c("rook", "bishop", "queen"),
  regional.weight = c("uniform", "frequency", "range"),
  abundance = FALSE,
  scale01 = TRUE
)
```

## Arguments

- pd:

  An object of class `vilma.pd` returned by
  [`rao.calc()`](https://oleon12.github.io/vilma/reference/rao.calc.md),
  containing the observed per-cell RaoQ table (`$pd.table`).

- tree:

  A rooted phylogeny of class `phylo` with branch lengths.

- dist:

  A `vilma.dist` object (e.g., from `points.to.raster()`) that includes
  `$distribution` with columns `Cell`, `Sp`, and `$grid` used when
  rasterizing per-cell results.

- iterations:

  Integer; number of randomizations (default `999`).

- method:

  Character; whether to compute null statistics at the `"global"` level
  (summing RaoQ across cells) or per `"cell"` (default behavior resolves
  to `"global"` if multiple are supplied).

- sampling:

  Character; null model type:

  - `"taxa.label"` — Shuffle tip labels on the tree; community matrix
    fixed.

  - `"range"` — Swap occurrences in a presence–absence matrix (preserves
    row/column sums).

  - `"neighbor"` — Local swaps between adjacent cells on `dist$grid`.

  - `"regional"` — Sample species for each cell from the regional pool
    (see `regional.weight`).

  If multiple options are passed, defaults to `"taxa.label"`.

- n.directions:

  Neighborhood definition for the `"neighbor"` model: `"rook"`,
  `"bishop"`, or `"queen"` (default resolves to `"queen"`).

- regional.weight:

  Weighting of the regional pool for the `"regional"` model: `"uniform"`
  (equal), `"frequency"` (by number of records), or `"range"` (by number
  of occupied cells). Defaults to `"uniform"`.

- abundance:

  Logical; if `FALSE` (default) each Cell×Sp table is binarized prior to
  computing proportions; if `TRUE` the counts from `table(Cell, Sp)` are
  used as abundances.

- scale01:

  Logical; if `TRUE` (default) patristic distances are divided by their
  maximum off-diagonal so \\D \in \[0,1\]\\ (diagonal set to 0).

## Value

An object of class `vilma.null`:

- `pd.obs`:

  Observed total RaoQ (global mode) or `CellValues$PD` (cell mode).

- `null.pd`:

  Vector of null totals (global mode) or omitted (cell mode).

- `SES`:

  Standardized effect size (scalar in global mode; per-cell in cell
  mode).

- `Pvalue`:

  Monte Carlo p-value(s).

- `CellValues`:

  Data frame with per-cell results (cell mode).

- `Raster`:

  Raster of SES (cell mode).

- `Iterations`:

  Number of randomizations used.

- `Iter.table`:

  Matrix of null RaoQ values (rows = cells, cols = iterations).

- `Method`:

  The string `"global"` or `"cell"`.

## Details

For each iteration, a randomized dataset (or tree, for `"taxa.label"`)
is produced, and
[`rao.calc()`](https://oleon12.github.io/vilma/reference/rao.calc.md) is
evaluated to obtain a null RaoQ per cell. The function returns either:

- **Global**: observed sum of RaoQ across cells, the vector of null
  sums, and `SES = (obs - mean(null)) / sd(null)` with a one-sided Monte
  Carlo p-value \\( \\(null \ge obs) + 1 ) / (N + 1)\\.

- **Cell**: per-cell observed RaoQ, null mean and SD, SES per cell, and
  per-cell p-values computed the same way, plus a raster of SES.

Species are restricted to the intersection between `tree$tip.label` and
`dist$distribution$Sp`. For `"neighbor"`, adjacency is derived from
`dist$grid` and the chosen `n.directions`. For `"regional"`, the number
of species sampled in each cell equals its observed richness
(`pd$pd.table$SR`), with probabilities defined by `regional.weight`.

## Interpretation

Negative SES indicates **phylogenetic clustering** (lower RaoQ than
expected), positive SES indicates **overdispersion** (higher RaoQ than
expected), and SES near 0 suggests random assembly under the chosen null
model.

## References

Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified
approach. *Theoretical Population Biology*, 21, 24–43.  
Botta-Dukát, Z. (2005). Rao’s quadratic entropy as a measure of
functional diversity. *Community Ecology*, 6, 283–290.  
Hardy, O. J. (2008). Testing the spatial phylogenetic structure of local
communities: statistical performances of different null models and test
statistics on a locally neutral community. *Journal of Ecology*, 96(5),
914–926.

## See also

[`rao.calc`](https://oleon12.github.io/vilma/reference/rao.calc.md) for
observed α-diversity (RaoQ);
[`rao.beta`](https://oleon12.github.io/vilma/reference/rao.beta.md) for
Rao β dissimilarity;
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)
for building `vilma.dist`.

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{

tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)

rao <- rao.calc(tree, dist)

out.g <- rao.calc.null(rao, tree = tree, dist = dist,
                       iterations = 999, method = "global",
                       sampling = "taxa.label")
print(out.g)
plot(out.g)
view.vilma(out.g)

# Per-cell SES under regional pool (frequency weights)
out.c <- rao.calc.null(pd = pd_obs, tree = tree, dist = dist,
                       iterations = 499, method = "cell",
                       sampling = "regional", regional.weight = "frequency")
} # }
```
