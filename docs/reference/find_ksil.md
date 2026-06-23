# Select optimal K using silhouette width

Computes hierarchical clustering on a phylogenetic dissimilarity matrix
and evaluates the mean silhouette width for a range of K values. The
function returns the K value that maximizes the average silhouette
width.

## Usage

``` r
find_ksil(pd.dis, method = c("ward.D2", "average", "complete"))
```

## Arguments

- pd.dis:

  A `dist` object representing pairwise dissimilarities among cells.

- method:

  Character string indicating the agglomeration method to be used in
  `hclust`. Options are `"ward.D2"`, `"average"`, or `"complete"`.
  Default is `"ward.D2"`.

## Value

A list containing:

- `best_k` The K value with the highest mean silhouette width.

- `sil_values` A data.frame with tested K values and their corresponding
  average silhouette widths.

## Examples

``` r
if (FALSE) { # \dontrun{
d <- dist(matrix(rnorm(100), 20, 5))
res <- find_ksil(d)
res$best_k
} # }
```
