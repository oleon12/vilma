#' Compute children list once (parent -> vector of children)
#' @keywords internal
#' @noRd

.children_list <- function(tree){
  split(tree$edge[,2], tree$edge[,1])
}

# Get descendant tips of an internal node using an iterative stack
.tips_from_node <- function(children, ntip, node){
  out <- integer(0)
  stack <- node
  while(length(stack)){
    v <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (v <= ntip) {
      out <- c(out, v)           # it's a tip index
    } else {
      ch <- children[[as.character(v)]]
      if (length(ch)) stack <- c(stack, ch)
    }
  }
  out
}

# Return a,b,c branch-length sums for two communities
branch.partition <- function(tree, commA, commB){
  ntip <- length(tree$tip.label)
  children <- .children_list(tree)

  # descendant tip names for each edge's child node
  desc_tips <- lapply(tree$edge[,2], function(node){
    idx <- if (node <= ntip) node else .tips_from_node(children, ntip, node)
    tree$tip.label[idx]
  })

  presA <- vapply(desc_tips, function(t) any(t %in% commA), TRUE)
  presB <- vapply(desc_tips, function(t) any(t %in% commB), TRUE)

  L <- tree$edge.length
  a <- sum(L[presA &  presB])
  b <- sum(L[presA & !presB])
  c <- sum(L[!presA &  presB])

  c(blShare = a, blA = b, blB = c)
}
