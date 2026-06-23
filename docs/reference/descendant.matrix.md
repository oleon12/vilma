# Build Descendant Matrix (edges × tips) for a Phylogenetic Tree

Construct a binary matrix `D` indicating which tips (species) descend
from the *child node* of each edge in a rooted phylogenetic tree. This
representation is useful for edge-based metrics (e.g., UniFrac) and
lineage-weighted beta-diversity decompositions.

## Usage

``` r
descendant.matrix(tree)
```

## Arguments

- tree:

  A rooted phylogenetic tree of class `phylo`. Branch lengths are not
  required for this function.

## Value

An integer matrix with:

- **Rows** corresponding to edges (`nrow(tree$edge)`), in the same order
  as `tree$edge`. Row names are `"parent-child"` node IDs (e.g.,
  `"14-15"`).

- **Columns** corresponding to tips (`tree$tip.label`).

- **Entries** `D[e, t] = 1` if tip `t` descends from the child node of
  edge `e`; otherwise `0`.

Properties:

- Terminal (tip) edges have exactly one `1`.

- Internal edges have as many `1`s as there are tips in the clade
  subtended by that edge.

- The output is deterministic and does not depend on edge reordering.

## Details

The algorithm computes, for every node, the set of descendant tips by
aggregating children → parent (a tips-to-root traversal). For each edge,
the descendant set of its child node defines which tip columns receive
1s.

## Note

The tree must be *rooted*. If your tree is unrooted, root it (e.g., with
[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html)) before calling
this function.

\#' @author Omar Daniel Leon-Alvarado
<https://leon-alvarado.weebly.com/> J. Angel Soto-Centeno
<https://www.mormoops.com/>

## See also

[`ape::is.rooted`](https://rdrr.io/pkg/ape/man/root.html),
[`ape::root`](https://rdrr.io/pkg/ape/man/root.html),
[`ape::read.tree`](https://rdrr.io/pkg/ape/man/read.tree.html)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ape)

tree <- data(Tree)

D <- descendant.matrix(tre)
dim(D)             # n_edge × n_tip
colnames(D)        # should match tr$tip.label
table(rowSums(D))  # 1 for terminals; >1 for internal edges
} # }

```
