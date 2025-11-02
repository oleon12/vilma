#' Find the Root Node of a Rooted Phylogeny
#' @description
#' This function identifies the numerical identifier of the root node in a
#' phylogenetic tree of class \code{phylo} from the \code{ape} package.
#' It is a helper function for navigating the structure of \code{phylo} objects.
#'
#' @param tree A phylogenetic tree. Must be an object of class \code{phylo}
#'   from the \code{ape} package and must be rooted. The function will throw an
#'   error if an unrooted tree is provided.
#'
#' @return An integer value corresponding to the node number of the root.
#'   For a standard \code{phylo} object, this is typically the highest-numbered
#'   internal node (e.g., for a tree with 5 tips, the root is node 6).
#'
#' @examples
#' # Create a simple rooted tree
#' library(ape)
#' tree <- read.tree(text = "((A:1, B:1):1, C:2);")
#' find.root.node(tree) # Should return the root node number (e.g., 5)
#'
#' @seealso The function \code{link[ape]{is.rooted}}
#'
#' @export

find.root.node <- function(tree) {
  
  if (!is.rooted(tree)) {
    
    stop("The tree is not rooted.")
  
  }
  
  # The root node is the one that appears in the parent column (edge[,1])
  # but never appears in the child column (edge[,2])
  
  all_nodes <- unique(tree$edge[, 1]) # All parent nodes
  children <- unique(tree$edge[, 2])  # All child nodes
  
  # Find the parent node that is not a child of any other node
  root_node <- all_nodes[!all_nodes %in% children]
  
  # It should be a single node. If not, something is wrong.
  if (length(root_node) != 1) {
    
    stop("Could not find a unique root node. Is the tree valid?")
  
  }
  
  return(root_node)

}
