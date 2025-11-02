#' Calculate Faith's PD for Individual Taxa
#' @description
#' Computes the phylogenetic diversity (PD) associated with one or more taxa in a phylogenetic tree.
#' Two methods are available:
#' \itemize{
#'   \item \strong{root} – PD from each taxon to the root of the tree.
#'   \item \strong{node} – PD as the branch length leading to each taxon.
#' }
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo}.
#' @param sp A character vector of tip labels to calculate PD for. Defaults to all tips in the tree.
#' @param method Character, either \code{"root"} or \code{"node"} indicating the PD calculation method. Defaults to \code{"root"}.
#'
#' @details
#' This function calculates Faith's PD for individual taxa. 
#' The \code{"root"} method sums branch lengths from the taxon to the root of the tree, 
#' while the \code{"node"} method returns only the branch leading directly to the taxon.
#' 
#' @return A named numeric vector of PD values for the specified taxa.
#'
#' @author Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#' 
#' @export

pd.taxon <- function(tree = NULL, sp = NULL, method = c("root","node")){
  
  # Only root trees are acepted Control 1
  if (!is.rooted(tree)) {
    stop("The tree is not rooted.")
  }
  
  #If no species is specify, then calculate all in the tree
  if(is.null(sp)==TRUE){
    sp <- tree$tip.label
  }
  
  #If no method is selected, then use root by default 
  if(length(method)==2){
    method <- "root"
    warning("method root by default")
  }
  
  ####################################

  spL <- 1:length(tree$tip.label)
  
  if(method == "root"){
    
    pd.out <- c()
    
    for(i in 1:length(sp)){
      
      #Get numeric position of each taxon
      spTMP <- spL[which(tree$tip.label%in%sp[i])]
      
      #Set stop and root nodes
      nodeStop <- spTMP
      nodeRoot <- find.root.node(tree)
      
      #Set the branch length object
      brl <- c()
      
      #This whill only stop once the nodeStop is the same as the node root
      while(nodeStop!=nodeRoot){
        
        #On child column, which position is the tip
        pos1 <- which(tree$edge[,2]%in%nodeStop)
        #Now, get the branch length of that position
        brl <- c(brl, tree$edge.length[pos1])
        #No move to the next node, so parental column will be
        #the new position, and restart the search at pos1
        nodeStop <- tree$edge[pos1,1]
        
      }
      
      #Once all branch length from one taxon are found, then sum them (pd)
      pd.out <- c(pd.out, sum(brl))
      
    }
    
    
  }
  
  if(method == "node"){
    
    pd.out <- c()
    
    for(i in 1:length(sp)){
      
      #Get numeric position of each taxon
      spTMP <- spL[which(tree$tip.label%in%sp[i])]
      
      #Set the branch length object
      brl <- c()
      
      #On child column, which position is the tip
      pos1 <- which(tree$edge[,2]%in%spTMP)
      
      #This method only take into account the branch length
      #of each taxon, thus, we save each single value as its pd
      pd.out <- c(pd.out, tree$edge.length[pos1])
      
    }
    
  }
  
  names(pd.out) <- sp
  
  return(pd.out)
}

