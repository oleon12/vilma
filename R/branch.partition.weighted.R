#' Precompute once per tree (outside loops) if you like:
#' @keywords internal
#' @noRd

.children_list <- function(tree) split(tree$edge[,2], tree$edge[,1])
.tips_from_node <- function(children, ntip, node){
  out <- integer(0); stack <- node
  while(length(stack)){
    v <- stack[[length(stack)]]; stack <- stack[-length(stack)]
    if (v <= ntip) out <- c(out, v) else {
      ch <- children[[as.character(v)]]; if (length(ch)) stack <- c(stack, ch)
    }
  }
  out
}
.desc_tips_by_edge <- function(tree){
  ntip <- length(tree$tip.label); children <- .children_list(tree)
  lapply(tree$edge[,2], function(node){
    idx <- if (node <= ntip) node else .tips_from_node(children, ntip, node)
    idx
  })
}

# Weighted branch partition: commA_abund / commB_abund are *named numeric* (names = tip labels)
branch.partition.weighted <- function(tree, desc_tips, commA_abund, commB_abund, normalize = TRUE){
  # Build tip-level abundance vectors aligned to tree$tip.label
  tips <- tree$tip.label
  a_vec <- rep(0, length(tips)); names(a_vec) <- tips
  b_vec <- rep(0, length(tips)); names(b_vec) <- tips

  if (length(commA_abund)) a_vec[names(commA_abund)] <- commA_abund
  if (length(commB_abund)) b_vec[names(commB_abund)] <- commB_abund

  if (normalize) {
    sa <- sum(a_vec)
    sb <- sum(b_vec)
    if (sa > 0) a_vec <- a_vec / sa
    if (sb > 0) b_vec <- b_vec / sb
  }

  # Edge-level abundances (sum of descendant tips)
  pA <- vapply(desc_tips, function(idx) sum(a_vec[idx]), numeric(1))
  pB <- vapply(desc_tips, function(idx) sum(b_vec[idx]), numeric(1))

  L <- tree$edge.length
  a <- sum(L * pmin(pA, pB))
  b <- sum(L * pmax(pA - pB, 0))
  c <- sum(L * pmax(pB - pA, 0))

  c(blShare = a, blA = b, blB = c)
}
