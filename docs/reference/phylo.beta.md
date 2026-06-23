# Phylogenetic Beta Diversity Partitioning (β_sor, β_sim, β_sne)

Computes pairwise **phylogenetic** beta diversity between communities
(cells) by partitioning total dissimilarity into **turnover** and
**nestedness** components, using branch-length–based analogs of
Baselga's decomposition.

Two modes are available:

- `method = "unweighted"` – presence/absence on the tree.

- `method = "weighted"` – relative-abundance weighting along edges.

Optional `normalize = TRUE` (default) rescales per-community abundances
to sum to 1 before edge aggregation, so indices lie in \\\[0,1\]\\ and
are comparable across communities with different total abundances.

## Usage

``` r
phylo.beta(tree, dist, method = c("unweighted", "weighted"), normalize = TRUE)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths.

- dist:

  An object of class `vilma.dist` containing species-by-cell occurrences
  (columns typically: `Sp`, `Lon`, `Lat`, `Cell`) and a reference grid.

- method:

  Character, either `"unweighted"` (default) or `"weighted"`. See
  *Mathematics* for definitions.

- normalize:

  Logical; if `TRUE` (default) and `method = "weighted"`, converts
  within-cell abundances to relative abundances (sum to 1) prior to edge
  aggregation. Ignored for `method = "unweighted"`.

## Value

An object of class `vilma.beta`:

- `total.dissimilarity` – `dist` object of \\\beta\_{\mathrm{sor}}\\.

- `turnover` – `dist` object of \\\beta\_{\mathrm{sim}}\\.

- `nestedness` – `dist` object of \\\beta\_{\mathrm{sne}}\\.

- `rasters` – named list with:

  - `mean.total`, `mean.turnover`, `mean.nestedness`.

  - PCoA rasters: `total.pcoa.1`, `total.pcoa.2`, `turnover.pcoa.1`,
    `turnover.pcoa.2`, `nestedness.pcoa.1`, `nestedness.pcoa.2`.

  - NMDS rasters: `total.ndms.1`, `total.ndms.2`, `turnover.ndms.1`,
    `turnover.ndms.2`, `nestedness.ndms.1`, `nestedness.ndms.2`.

- `stats` – data frame summarizing eigenvalues (PCoA) and stress (NMDS).

- `calculation.method` – resolved `method` and whether normalization was
  applied.

- `algorithm` – the string `"PhyloBeta"`.

## Details

Only species present in both `tree` and `dist` are used; informative
messages report mismatches. Pairwise dissimilarities are computed for
all cell combinations. Ordinations (PCoA with Cailliez correction; NMDS)
are run on each dissimilarity to obtain spatial axes, and per-cell means
are computed by averaging off-diagonal dissimilarities.

## Mathematics

Let \\E\\ be the set of tree edges with lengths \\L_e\\. For two
communities (cells) \\A\\ and \\B\\, define for each edge \\e\\ the
fraction of the community that lies *below* that edge:

- *Unweighted (presence/absence):* \\p_A(e) = \mathbf{1}\\\text{any
  descendant tip of } e \text{ in } A\\\\, and similarly \\p_B(e)\\.

- *Weighted (relative abundance):* \\p_A(e) = \sum\_{i \in
  \text{tips}(e)} \tilde{a}\_i\\, where \\\tilde{a}\_i = a_i / \sum_j
  a_j\\ are per-tip relative abundances in community \\A\\ (and
  analogously \\p_B(e)\\ for \\B\\).

From these, define branch-length partitions: \$\$a = \sum\_{e \in E} L_e
\\ \min\\p_A(e), p_B(e)\\,\$\$ \$\$b = \sum\_{e \in E} L_e \\
\max\\p_A(e) - p_B(e),\\ 0\\, \qquad c = \sum\_{e \in E} L_e \\
\max\\p_B(e) - p_A(e),\\ 0\\.\$\$ The phylogenetic beta components are
then \$\$\beta\_{\mathrm{sor}} = \frac{b + c}{\\2a + b + c\\}, \qquad
\beta\_{\mathrm{sim}} = \frac{\min(b,c)}{\\a + \min(b,c)\\}, \qquad
\beta\_{\mathrm{sne}} = \beta\_{\mathrm{sor}} -
\beta\_{\mathrm{sim}}.\$\$ These reduce to Baselga's (2010) PA formulas
when \\p_A(e), p_B(e) \in \\0,1\\\\. By construction \\0 \le
\beta\_{\mathrm{sim}}, \beta\_{\mathrm{sne}}, \beta\_{\mathrm{sor}} \le
1\\ and \\\beta\_{\mathrm{sor}} = \beta\_{\mathrm{sim}} +
\beta\_{\mathrm{sne}}\\.

## Implementation Notes

- Presence/absence mode is obtained by setting \\p_A(e), p_B(e)\\ to
  edge-level indicators of occupancy.

- In weighted mode, setting `normalize = TRUE` yields relative-abundance
  fractions (sum to 1 per cell) before edge aggregation; `FALSE` uses
  raw counts and can emphasize differences in total sampling intensity.

- The distance matrices are symmetric with zero diagonals; internal
  checks confirm bounds and \\\beta\_{\mathrm{sor}} =
  \beta\_{\mathrm{sim}} + \beta\_{\mathrm{sne}}\\.

- A progress bar can be used during pairwise loops; when present, it
  should be closed with `on.exit(close(pb), add = TRUE)`.

## References

Baselga, A. (2010). Partitioning the turnover and nestedness components
of beta diversity. *Global Ecology and Biogeography*, 19, 134–143.
Leprieur, F., Albouy, C., de Bortoli, J., Cowman, P.F., Bellwood, D.R.,
& Mouillot, D. (2012). Quantifying phylogenetic beta diversity:
distinguishing turnover of lineages from PD gradients. *PLoS ONE*, 7(8),
e42760. Nipperess, D.A., & Matsen IV, F.A. (2013). The mean and variance
of phylogenetic diversity under rarefaction. *Methods in Ecology and
Evolution*, 4, 566–572.

## See also

[`unifrac.calc`](https://oleon12.github.io/vilma/reference/unifrac.calc.md),
[`phylosor.calc`](https://oleon12.github.io/vilma/reference/phylosor.calc.md),
[`faith.pd`](https://oleon12.github.io/vilma/reference/faith.pd.md),
[`points_to_raster`](https://oleon12.github.io/vilma/reference/points_to_raster.md)

## Author

Omar Daniel Leon-Alvarado <https://leon-alvarado.weebly.com/> J. Angel
Soto-Centeno <https://www.mormoops.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- examplae_tree()
dist <- example_dist()
dist <- points_to_raster(dist)

# Presence/absence (Baselga-style) phylogenetic β
pb_unw <- phylo.beta(tree, dist, method = "unweighted")
print(pb_unw)
plot(pb_unw)
view.vilma(pb_unw)

# Abundance-weighted phylogenetic β (relative abundances by default)
pb_w <- phylo.beta(tree, dist, method = "weighted", normalize = TRUE)
} # }
```
