# Find the Root Node of a Rooted Phylogeny

This function identifies the numerical identifier of the root node in a
phylogenetic tree of class `phylo` from the `ape` package. It is a
helper function for navigating the structure of `phylo` objects.

## Usage

``` r
find.root.node(tree)
```

## Arguments

- tree:

  A phylogenetic tree. Must be an object of class `phylo` from the `ape`
  package and must be rooted. The function will throw an error if an
  unrooted tree is provided.

## Value

An integer value corresponding to the node number of the root. For a
standard `phylo` object, this is typically the highest-numbered internal
node (e.g., for a tree with 5 tips, the root is node 6).

## See also

The function `link[ape]{is.rooted}`

## Examples

``` r
# Create a simple rooted tree
library(ape)
tree <- read.tree(text = "((A:1, B:1):1, C:2);")
find.root.node(tree) # Should return the root node number (e.g., 5)
#> [1] 4
```
