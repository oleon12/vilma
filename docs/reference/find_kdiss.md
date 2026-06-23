# Select optimal K using explained dissimilarity

Computes hierarchical clustering on a dissimilarity matrix and evaluates
explained dissimilarity (R2) across a range of K values. Explained
dissimilarity is computed as `R2 = 1 - (W/T)`, where `W` is the sum of
within-cluster pairwise dissimilarities and `T` is the sum of all
pairwise dissimilarities. The function selects an optimal K using either
a second-difference ("2diff") knee rule or a relative threshold ("thr")
based on the incremental gain in R2.

## Usage

``` r
find_kdiss(
  pd.dis,
  method = c("ward.D2", "average", "complete"),
  k_thr = c("2diff", "thr"),
  thr = 0.1
)
```

## Arguments

- pd.dis:

  A `dist` object representing pairwise dissimilarities among cells.

- method:

  Character string indicating the agglomeration method to be used in
  `hclust`. Options are `"ward.D2"`, `"average"`, or `"complete"`.
  Default is `"ward.D2"`.

- k_thr:

  Character string indicating the K-selection rule. Options are
  `"2diff"` (knee via second difference) or `"thr"` (relative
  threshold). Default is `"2diff"`.

- thr:

  Numeric value used only when `k_thr = "thr"`. Represents the
  proportion of the maximum incremental gain in R2 used as a cutoff. The
  function selects the smallest K such that `dR2 < thr * max(dR2)`. For
  example, `thr = 0.1` selects the first K where the gain in explained
  dissimilarity is less than 10% of the maximum observed gain.

## Value

A list containing:

- `best_k` The selected K value.

- `R_vals` A data.frame with tested K values, `R2.vals`, and `dR2`.

## Examples

``` r
if (FALSE) { # \dontrun{
d <- dist(matrix(rnorm(100), 20, 5))
res <- find_kdiss(d, method = "ward.D2", k_thr = "2diff")
res$best_k
head(res$R_vals)
} # }
```
