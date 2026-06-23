# Rao Beta (phylogenetic) dissimilarity index

Computes **between-community** phylogenetic dissimilarity (Rao β) among
spatial cells using a rooted phylogeny and a `vilma.dist` object. Rao β
is the expected phylogenetic dissimilarity between two individuals drawn
from different communities: \$\$\Delta\_{\mathrm{Rao}}(i,j) =
-\tfrac{1}{2}\left(Q_i + Q_j - 2\\p_i^\top D\\p_j\right),\$\$ where \\Q
= p^\top D p\\ is Rao’s α within a community, \\p\\ are within-cell
relative abundances (presence–absence if `abundance = FALSE`), and \\D\\
is the patristic distance matrix among species.

## Usage

``` r
rao.beta(tree, dist, abundance = FALSE, scale01 = TRUE)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths.

- dist:

  An object of class `vilma.dist` (e.g., from `points.to.raster()`),
  containing at least columns `Cell` and `Sp` in `$distribution`, and a
  raster `$grid` for mapping outputs.

- abundance:

  Logical. If `FALSE` (default), the Cell×Sp table is binarized
  (presence–absence) before proportions are computed. If `TRUE`, raw
  Cell×Sp counts from `table(Cell, Sp)` define relative abundances.

- scale01:

  Logical. If `TRUE` (default), \\D\\ is divided by its maximum
  off-diagonal so \\D \in \[0,1\]\\ (diagonal set to 0). Under this
  scaling, Rao β lies in \\\[0,1\]\\ up to numerical tolerance.

## Value

An object of class `vilma.beta`:

- `distribution` – Original `dist$distribution`.

- `Rao.Beta` – Pairwise Rao β dissimilarity as a
  [`stats::dist`](https://rdrr.io/r/stats/dist.html) object.

- `rasters` – List with: `mean.dissimilarity`, `pcoa.1`, `pcoa.2`,
  `ndms.1`, `ndms.2`.

- `pcoa.eig` – Relative corrected eigenvalues for axes 1–2.

- `ndms.stress` – NDMS stress value.

- `calculation.method` – Summary string of `abundance` and `scale01`.

- `algorithm` – The string `"RaoBeta"`.

## Details

Species are intersected between `tree$tip.label` and
`dist$distribution$Sp`; only the common set is used. To ensure correct
behavior under taxa-label null models, the patristic matrix \\D\\ is
**reindexed to match** the community matrix column order before matrix
multiplications. Empty cells (no species after optional binarization)
have their corresponding rows/columns of the output dissimilarity set to
`NA`.

In addition to the full pairwise Rao β matrix, the function returns:

- the mean Rao β per cell (raster),

- a PCoA of the dissimilarity (first two axes rasterized; Cailliez
  correction),

- a 2D NDMS (first two axes rasterized) and the stress value.

## Interpretation

With `scale01 = TRUE`, Rao β ranges from `0` (identical phylogenetic
composition) to `1` (maximal turnover). A similarity matrix can be
obtained as `1 - as.matrix(Rao.Beta)` if needed.

## References

Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified
approach. *Theoretical Population Biology*, 21, 24–43.  
Botta-Dukát, Z. (2005). Rao’s quadratic entropy as a measure of
functional diversity. *Community Ecology*, 6, 283–290.

## See also

[`rao.calc`](https://oleon12.github.io/vilma/reference/rao.calc.md) (Rao
α),
[`rao.calc.null`](https://oleon12.github.io/vilma/reference/rao.calc.null.md)
(null models),
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
out <- rao.beta(tree = tree, dist = dist, abundance = FALSE, scale01 = TRUE)
print(out)
plot(out)
view.vilma(out)
} # }
```
