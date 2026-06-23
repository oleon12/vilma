# UniFrac Phylogenetic Dissimilarity between communities (cells)

Computes pairwise UniFrac distances among grid cells and returns: (i) a
dissimilarity matrix, and (ii) convenient rasters summarizing mean
dissimilarity per cell as well as ordination (PCoA and NMDS) axes mapped
to space.

The computation follows the classic UniFrac definitions: when
`method = "unweighted"`, distances depend only on presence/absence; when
`method = "weighted"`, distances are normalized and depend on the
relative abundances of taxa in each cell.

## Usage

``` r
unifrac.calc(tree, dist, method = c("unweighted", "weighted"))
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo` with branch lengths.

- dist:

  An object of class `vilma.dist` (e.g., from `points.to.raster()`),
  containing species occurrences per grid cell and a template raster.

- method:

  Character; UniFrac variant to compute:

  - `"unweighted"` – presence/absence UniFrac (Default).

  - `"weighted"` – normalized weighted UniFrac using relative
    abundances.

## Value

An object of class `vilma.beta`, a list with:

- `distribution` – Original species distribution (from `dist`).

- `UniFrac` – UniFrac dissimilarity as a `dist` object (labels are cell
  IDs).

- `rasters` – A named list of `SpatRaster` layers:

  - `mean.unifrac` – Raster of mean UniFrac per cell (diagonal
    excluded).

  - `pcoa.1`, `pcoa.2` – PCoA axes 1 and 2 mapped to cells (from
    dissimilarity).

  - `ndms.1`, `ndms.2` – NMDS axes 1 and 2 mapped to cells (from
    dissimilarity).

- `pcoa.eig` – PCoA eigen values of the two first axes.

- `ndms.stress` – NMDS stress value.

- `calculation.method` – The resolved `method` used.

- `algorithm` – The string `"UniFrac"`.

## Details

Let \\E\\ be the set of tree edges with lengths \\w_e\\. For two
communities (cells) \\A\\ and \\B\\, define \\U_e(A)\\ and \\U_e(B)\\ as
indicators that edge \\e\\ has at least one descendant tip present in
community \\A\\ (or \\B\\), respectively.

**Unweighted UniFrac** (presence/absence) is:
\$\$\mathrm{UniFrac}\_{\mathrm{unw}}(A,B) = \frac{\sum\_{e \in E} w_e
\cdot \mathbf{1}\\ U_e(A) + U_e(B) = 1 \\} {\sum\_{e \in E} w_e \cdot
\mathbf{1}\\ U_e(A) + U_e(B) \ge 1 \\}.\$\$

**Weighted UniFrac** (normalized; relative abundances) defines the
fraction of abundance below each edge as \\D_A(e)\\ and \\D_B(e)\\ (each
in \\\[0,1\]\\), and uses: \$\$\mathrm{UniFrac}\_{\mathrm{w}}(A,B) =
\frac{\sum\_{e \in E} w_e \\ \|D_A(e) - D_B(e)\|} {\sum\_{e \in E} w_e
\\ \[D_A(e) + D_B(e)\]}.\$\$

Only species shared between `tree` and `dist` are used; informative
messages are printed for mismatches. The returned distance matrix has
diagonal 0.

After computing the matrix, the function summarizes each cell by its
mean dissimilarity to all other cells and performs ordinations on the
dissimilarity matrix using PCoA (with Cailliez correction) and NMDS,
returning the first two axes rasterized to the study grid.

## Notes

- The tree must be rooted and have branch lengths.

- Only species present in both `tree` and `dist` are used.

- Ordinations are performed on the dissimilarity matrix. PCoA uses
  `ape::pcoa(..., correction = "cailliez")` to handle non-Euclidean
  cases.

- `method = "weighted"` implements the *normalized* weighted UniFrac
  (values in \\\[0,1\]\\), based on relative abundances.

## References

Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method
for comparing microbial communities. *Applied and Environmental
Microbiology*, 71(12), 8228–8235.

Lozupone, C., Hamady, M., Kelley, S.T., & Knight, R. (2007).
Quantitative and qualitative beta diversity measures lead to different
insights into factors that structure microbial communities. *Applied and
Environmental Microbiology*, 73(5), 1576–1585.

Chen, J., Bittinger, K., Charlson, E.S., et al. (2012). Associating
microbiome composition with environmental covariates using generalized
UniFrac distances. *Bioinformatics*, 28(16), 2106–2113.

## See also

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

# Unweighted UniFrac (presence/absence)
uf_unw <- unifrac.calc(tree, dist, method = "unweighted")
print(uf_unw)
plot(uf_unw)
view.vilma(uf_unw)

# Weighted (normalized) UniFrac (relative abundances)
uf_w <- unifrac.calc(tree, dist, method = "weighted")
} # }

```
