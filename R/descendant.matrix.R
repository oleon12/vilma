#' Build Descendant Matrix (edges × tips) for a Phylogenetic Tree
#'
#' Construct a binary matrix \code{D} indicating which tips (species) descend
#' from the \emph{child node} of each edge in a rooted phylogenetic tree.
#' This representation is useful for edge-based metrics (e.g., UniFrac) and
#' lineage-weighted beta-diversity decompositions.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo}.
#'   Branch lengths are not required for this function.
#'
#' @details
#' The algorithm computes, for every node, the set of descendant tips by
#' aggregating children → parent (a tips-to-root traversal). For each edge,
#' the descendant set of its child node defines which tip columns receive 1s.
#' 
#' @return An integer matrix with:
#' \itemize{
#'   \item \strong{Rows} corresponding to edges (\code{nrow(tree$edge)}), in the
#'         same order as \code{tree$edge}. Row names are \code{"parent-child"}
#'         node IDs (e.g., \code{"14-15"}).
#'   \item \strong{Columns} corresponding to tips (\code{tree$tip.label}).
#'   \item \strong{Entries} \code{D[e, t] = 1} if tip \code{t} descends from the
#'         child node of edge \code{e}; otherwise \code{0}.
#' }
#'
#'
#' Properties:
#' \itemize{
#'   \item Terminal (tip) edges have exactly one \code{1}.
#'   \item Internal edges have as many \code{1}s as there are tips in the clade
#'         subtended by that edge.
#'   \item The output is deterministic and does not depend on edge reordering.
#' }
#'
#' @note The tree must be \emph{rooted}. If your tree is unrooted, root it
#' (e.g., with \code{ape::root()}) before calling this function.
#'
#' #' @author 
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com/}
#'
#' @seealso \code{ape::is.rooted}, \code{ape::root}, \code{ape::read.tree}
#'
#' @seealso \code{ape::is.rooted}, \code{ape::root}, \code{ape::read.tree}
#' 
#'@examples
#' \dontrun{
#' library(ape)
#'
#' tree <- data(Tree)
#'
#' D <- descendant.matrix(tre)
#' dim(D)             # n_edge × n_tip
#' colnames(D)        # should match tr$tip.label
#' table(rowSums(D))  # 1 for terminals; >1 for internal edges
#' }
#'
#'
#' @export

descendant.matrix <- function(tree) {
  if (!inherits(tree, "phylo")) stop("Input must be a 'phylo' object.")
  if (!is.rooted(tree)) stop("Tree must be rooted.")
  
  ntip  <- length(tree$tip.label)
  nedge <- nrow(tree$edge)
  
  # Initialize descendant list: each tip is its own descendant
  des <- vector("list", ntip + tree$Nnode)
  for (i in seq_len(ntip)){
    des[[i]] <- i
  }
  # Build list of children for every internal node
  children <- split(tree$edge[, 2], tree$edge[, 1])
  
  # Traverse nodes from tips upward (postorder)
  # So every child’s descendants are computed before the parent
  post <- rev(unique(tree$edge[, 1]))  # simple topological order works fine
  
  for (node in post) {
    if (node > ntip) {
      child_nodes <- children[[as.character(node)]]
      if (!is.null(child_nodes)) {
        des[[node]] <- sort(unique(unlist(des[child_nodes], use.names = FALSE)))
      }
    }
  }
  
  # Build D: rows = edges, columns = tips
  D <- matrix(0L, nrow = nedge, ncol = ntip,
              dimnames = list(apply(tree$edge, 1, paste, collapse = "-"),
                              tree$tip.label))
  
  for (e in seq_len(nedge)){
    
    child <- tree$edge[e, 2]
    tips  <- des[[child]]
    
    if (length(tips)){
      D[e, tips] <- 1L
    }
  }
  return(D)
}

