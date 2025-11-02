#' Null model for Faith's Phylogenetic Diversity (PD)
#' @description
#' This function implements null models for Faith's PD using different
#' randomization approaches. It allows testing significance of PD at either 
#' the global (across all cells) or cell level, with several sampling 
#' strategies: randomizing taxa labels, species ranges, or neighborhood swaps.
#'
#' @param pd An object of class \code{vilma.pd}, containing observed PD values.
#' @param tree A rooted phylogenetic tree of class \code{phylo}, with branch lengths.
#' @param dist An object of class \code{vilma.dist}, representing species distributions
#'   (see \code{points.to.raster()}).
#' @param iterations Integer. Number of randomizations to perform (default = 999).
#' @param method Character. Null model output method:
#'   \itemize{
#'     \item \code{"global"} – significance of total PD across all cells.
#'     \item \code{"cell"} – significance per cell.
#'   }
#' @param sampling Character. Randomization strategy:
#'   \itemize{
#'     \item \code{"taxa.label"} – permutes species labels on the tree.
#'     \item \code{"range"} – randomizes species ranges using swap algorithms.
#'     \item \code{"neighbor"} – swaps species occurrences between adjacent cells.
#'   }
#' @param n.directions Character. Neighborhood adjacency definition for the
#'   \code{"neighbor"} method. Options: \code{"rook"}, \code{"bishop"}, \code{"queen"}
#'   (default = \code{"queen"}).
#'
#' @param regional.weight Weighting of the regional pool for the \code{"regional"} model:
#'   \code{"uniform"} (equal), \code{"frequency"} (by number of records), or
#'   \code{"range"} (by number of occupied cells). Defaults to \code{"uniform"}.
#'
#' @details
#' The function generates a null distribution of PD values based on the chosen 
#' sampling strategy. For the \code{"global"} method, the observed PD is compared 
#' against the distribution of total PD values across iterations. For the 
#' \code{"cell"} method, standardized effect sizes (SES) and p-values are computed 
#' for each grid cell individually.
#'
#' @return An object of class \code{vilma.null}, with structure depending on \code{method}:
#' \itemize{
#'   \item For \code{global}: a list with observed PD, null distribution, SES, p-value,
#'         iteration table, and metadata.
#'   \item For \code{cell}: a list with per-cell values (PD, null mean, null SD, SES, p-value),
#'         a raster of SES values, iteration table, and metadata.
#' }
#'
#' @references
#' Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. 
#' \emph{Biological Conservation}, 61(1), 1–10. \doi{10.1016/0006-3207(92)91201-3}
#'
#' Gotelli, N. J., & Graves, G. R. (1996). \emph{Null Models in Ecology}. 
#' Smithsonian Institution Press.
#'
#' Gotelli, N. J. (2000). Null model analysis of species co‐occurrence patterns. 
#' \emph{Ecology}, 81(9), 2606–2621. \doi{10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2}
#'
#' Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002). 
#' Phylogenies and community ecology. \emph{Annual Review of Ecology and Systematics}, 
#' 33, 475–505. \doi{10.1146/annurev.ecolsys.33.010802.150448}
#'
#' @author 
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com/}
#'
#' @seealso \code{\link{faith.pd}}, \code{\link{points_to_raster}}, 
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' pd <- faith.pd(tree, dist, method = "cell")
#' null_model <- faith.pd.null(pd, tree, dist,
#'                             iterations = 999, method = "cell",
#'                             sampling = "range")
#' }
#'
#'
#' @export


faith.pd.null <- function(pd, tree, dist, iterations = 999,
                                   method = c("global","cell"),
                                   sampling = c("taxa.label","range","neighbor","regional"),
                                   n.directions = c("rook","bishop","queen"),
                                   regional.weight = c("uniform","frequency","range")
                                   ){
  
  ##############################################################
  #                        VERIFICATION                        #
  ##############################################################
  
  if (!inherits(pd, "vilma.pd")) {
    stop("Input 'pd' must be an object of class 'vilma.pd'.")
  }
  
  if(!inherits(tree, "phylo")){
    stop("Input 'tree' must be an objecto of class 'phylo'.")
    
  }
  
  if (is.null(tree$edge.length)) {
    stop("The phylogenetic tree must have branch lengths.")
  }
  
  if (!is.rooted(tree)) {
    stop("Faith's PD requires a rooted tree.")
  }
  
  if (!inherits(dist, "vilma.dist")) {
    stop("Input 'dist' must be an object of class 'vilma.dist'. See 'points_to_raster()' function.")
  }
  
  if(length(method) == 2){
    method <- "global"
  }
  
  if(length(sampling) != 1){
    sampling <- "taxa.label"
  }
  
  if(length(n.directions) != 1){
    n.directions <- "queen"
  }
  
  if(length(regional.weight) != 1){
    regional.weight <- "uniform"
  }
  
  ##############################################################
  #                      END VERIFICATION                      #
  ##############################################################
  
  #####################################################################
  
  #########################################################################
  
  if(sampling == "neighbor"){
    message( "Running null-model with: ", iterations, " iterations. 
	  \n", "Method: ", method, "\n", "Sampling: ", sampling, "\n", 
             "Neighbor direction: ", n.directions  )
  }else{
    if(sampling == "regional"){
      message("Running  null-model with: ", iterations, " iterations. \n",
              "Method: ", method, "\n", "Sampling: ", sampling, "\n",
              "Weight: ", regional.weight, "\n")
    }else{
      message( "Running null-model with: ", iterations, " iterations. \n", "Method: ", method, "\n", "Sampling: ", sampling, "\n"  )
    }
  }
  
  
  ##############################################################
  #                        doParallel                          #
  ##############################################################
  
  #########################################################################
  
  ###############################
  #         Taxa Method         #
  ###############################
  
  if(sampling == "taxa.label"){
    
    #Create the ouput table
    #rows = number of cells; cols = PD values per interaction
    samples <- matrix(NA, nrow = nrow(pd$pd.table), ncol = iterations)
    
    #Avoid overlapping and environment contamination
    
    dist1 <- dist
    
    #extract original tip.label
    taxa.label <- tree$tip.label
    
    #Take into account only species shared between the tree and dist
    taxa.label <- intersect(taxa.label, dist$distribution$Sp)
    #Find the original position of common spcies
    #Thus only those species will be resampled
    taxa.tree.pos <- which(tree$tip.label %in% taxa.label)
    
    
    set.seed(123)
    
    #Iterations
    
    iter.tree <- tree
    iter.tree.out <- list()
    
    message("\n","Iterating trees ..." ,"\n")
    
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    
    for(i in 1:iterations){
      iter.tree$tip.label[taxa.tree.pos] <- sample(taxa.label, length(taxa.label))
      iter.tree.out[[i]] <- iter.tree
      setTxtProgressBar(pb, i)
    }
    
    message("\n","Calculating PD ..." ,"\n")
    
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for(i in 1:iterations){
      
      #Calculate pd
      pd.iter <- NA
      pd.iter <- suppressMessages(faith.pd(iter.tree.out[[i]], dist1, method = pd$calculation.method))
      samples[,i] <- pd.iter$pd.table$PD
      
      
      setTxtProgressBar(pb, i)
      
    }
    
    #samples <- samples[,-1]
    colnames(samples) <- paste0("It.",as.character(1:ncol(samples)))
    
    
  }
  
  
  #########################################################################
  
  ###############################
  #        Range Method         #
  ###############################
  
  if(sampling  == "range"){
    
    #Create the ouput table
    #rows = number of cells; cols = PD values per interaction
    samples <- matrix(NA, nrow = nrow(pd$pd.table), ncol = iterations)
    
    #Avoid overlapping and environment contamination
    
    dist1 <- dist
    
    #Generate the presence/absence matrix for the swap
    
    #Using the original first
    dist0 <- dist$distribution[ , c(1,4)]
    
    #This mother fucker is the most important line
    #Doing this I guaranteed that only species shared between Dist and Tree
    #are the ones that are iterated, and avoid collapse any original cell
    dist01 <- dist0[which(dist0$Sp%in%intersect(dist0$Sp, tree$tip.label)), ]  
    
    #Turn it into a table
    dist01 <- table(dist01$Cell, dist01$Sp)
    dist01 <- (dist01 > 0) *1 #Only 1 & 0
    dist01 <- as.matrix(dist01)
    
    iter.dist <- list()
    
    message("\n","Swaping distributions ..." ,"\n")
    
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for(i in 1:iterations){
      iter.d <- swap.null(dist01)
      dist1$distribution <- return.vilma.dist(iter.d)
      iter.dist[[i]] <- dist1
      dist1$distribution <- NA
      setTxtProgressBar(pb, i)
      
    }
    
    message("\n","Calculating PD ..." ,"\n")
    
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for(i in 1:iterations){
      
      pd.iter <- NA
      pd.iter <- suppressMessages(faith.pd(tree, iter.dist[[i]], method = pd$calculation.method))
      
      samples[ , i] <- pd.iter$pd.table$PD
      
      setTxtProgressBar(pb, i)
      
    }
    
    #samples <- samples[,-1]
    colnames(samples) <- paste0("It.",as.character(1:ncol(samples)))
    
  }
  
  #########################################################################
  
  ###############################
  #       Neighbor Method       #
  ###############################
  
  if(sampling == "neighbor"){
    
    #Create the ouput table
    #rows = number of cells; cols = PD values per interaction
    samples <- matrix(NA, nrow = nrow(pd$pd.table), ncol = iterations)
    
    #This mother fucker is the most important line
    #Doing this I guaranteed that only species shared between Dist and Tree
    #are the ones that are iterated, and avoid collapse any original cell
    dist0 <- dist
    dist0$distribution <- dist0$distribution[which(dist0$distribution$Sp %in% intersect(dist0$distribution$Sp, tree$tip.label)), ] 
    
    dist0$distribution <- dist0$distribution[which(dist0$distribution$Cell %in% as.integer(pd$pd.table$Cell)), ]
    
    #Avoid overlapping and environment contamination
    dist1 <- dist0
    dist2 <- dist0
    
    #Iteration dist output
    iter.dist <- list()
    
    #Calculate the number of iterations for neighbor swap
    iter <- (ncell(dist1$grid)* (log(ncell(dist1$grid))+0.577))
    iter <- iter *2
    
    message("\n","Swaping neighbors ..." ,"\n")
    
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for(i in 1:iterations){
      #Ensure factors
      dist1$distribution$Sp <- as.factor(dist1$distribution$Sp)
      dist1$distribution$Cell <- as.integer(dist1$distribution$Cell)
      
      #Found adjacent cells
      n.cells <- adjacent(x = dist1$grid,
                                 cells = 1:ncell(dist1$grid),
                                 directions = n.directions,
                                 pairs = TRUE)
      #List, each cell with its respective adjacent cells
      n.cells <- split(n.cells[ ,2], n.cells[ ,1])
      
      #List, each cell with its respective spp.
      spp.cells <- split(dist1$distribution$Sp, dist1$distribution$Cell)
      
      for(j in 1:iter){
        
        # Pick one random cell
        cell.a <- sample(names(spp.cells),1)
        sp.a <- spp.cells[[cell.a]]
        if(length(sp.a) == 0) next
        
        # Pick one random neighbor
        neighs <- n.cells[[cell.a]]
        if(is.null(neighs) || length(neighs)==0) next
        cell.b <- as.character(sample(neighs,1))
        sp.b <- spp.cells[[cell.b]]
        if(length(sp.b) == 0) next
        
        # Pick one species from each cell
        sp.a.pick <- sample(sp.a, 1)
        sp.b.pick <- sample(sp.b, 1)
        
        # Only swap if species aren't already present
        if(!(sp.a.pick %in% sp.b) && !(sp.b.pick %in% sp.a)){
          
          # True swap: replace species in both cells
          spp.cells[[cell.a]][spp.cells[[cell.a]] == sp.a.pick][1] <- sp.b.pick
          spp.cells[[cell.b]][spp.cells[[cell.b]] == sp.b.pick][1] <- sp.a.pick
        }
      }
      
      Sp <- c()
      Cell <- c()
      
      for(x in 1:length(spp.cells)){
        Sp <- c(Sp, as.character(spp.cells[[x]]))
        Cell <- c(Cell, as.integer(rep(names(spp.cells)[x], length(spp.cells[[x]]))))
      }
      
      dist2$distribution <- NA
      dist2$distribution <- data.frame(Sp =Sp,
                                       Lon = rep(0, length(Sp)),
                                       Lat = rep(0, length(Sp)),
                                       Cell = Cell)
      iter.dist[[i]] <- dist2
      
      setTxtProgressBar(pb, i)
      
    }
    
    message("\n","Calculating PD ..." ,"\n")
    
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for(i in 1:iterations){
      
      pd.iter <- NA
      pd.iter <- suppressMessages(faith.pd(tree, iter.dist[[i]], method = pd$calculation.method))
      
      samples[ , i] <- pd.iter$pd.table$PD
      
      setTxtProgressBar(pb, i)
      
    }
    
    #samples <- samples[,-1]
    colnames(samples) <- paste0("It.",as.character(1:ncol(samples)))
    
  }
  
  #########################################################################
  
  ###############################
  #        Regional Pool        #
  ###############################
  
  if(sampling == "regional"){
    
    #Create the ouput table
    #rows = number of cells; cols = PD values per interaction
    samples <- matrix(NA, nrow = nrow(pd$pd.table), ncol = iterations)
    
    
    #This mother fucker is the most important line
    #Doing this I guaranteed that only species shared between Dist and Tree
    #are the ones that are iterated, and avoid collapse any original cell
    dist0 <- dist
    dist0$distribution <- dist0$distribution[which(dist0$distribution$Sp %in% intersect(dist0$distribution$Sp, tree$tip.label)), ] 
    
    dist0$distribution <- dist0$distribution[which(dist0$distribution$Cell %in% as.integer(pd$pd.table$Cell)), ]
    
    #To avoid environment contamination
    dist1 <- dist0
    
    #Preprocess species and weights
    distM <- dist0$distribution
    tree.sp <- tree$tip.label # Get species in phylo
    obs.sp <- unique(distM$Sp) # Get species in distribution
    common.sp <- intersect(tree.sp, obs.sp) #Use only unique species
    if (length(common.sp) == 0) stop("No overlapping species for regional null.")
    
    ## Uniform weight ##
    weight <- rep(1, length(common.sp))
    names(weight) <- common.sp
    
    ## Frequency weight ##
    # Number of cells each species occupies
    if(regional.weight == "frequency"){
      
      species.freq <- table(distM$Sp) # counts all rows = frequency
      weight[names(species.freq)] <- species.freq[names(species.freq)]
      
    }
    
    if(regional.weight == "range"){
      occ.tab <- table(distM$Sp, distM$Cell)
      species.cells <- rowSums(occ.tab > 0) # number of unique cells
      weight[names(species.cells)] <- species.cells[names(species.cells)]
      
    }
    
    # Normalize to probabilities
    prob_vec <- weight / sum(weight)
    
    # Observed cell richness
    pd_table <- pd$pd.table
    cell_ids <- pd_table$Cell
    richness_vec <- pd_table$SR
    
    #Iteration dist output
    iter.dist <- list()
    
    message("\n","Getting distributions ..." ,"\n")
    
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for(it in seq_len(iterations)){
      null_list <- vector("list", length = length(cell_ids))
      names(null_list) <- cell_ids
      
      for (i in seq_along(cell_ids)) {
        cid <- cell_ids[i]
        k <- richness_vec[i]
        if (k > 0) {
          # sample k species without replacement
          null_list[[i]] <- sample(common.sp, size = k, replace = FALSE, prob = prob_vec)
        } else {
          null_list[[i]] <- character(0)
        }
      }
      
      # convert null_list (species per cell) into a dist object
      # you need a helper function to build a vilma.dist-like object from null_list
      null_dist <- return.vilma.dist2(null_list)
      
      dist1$distribution <- null_dist
      
      iter.dist[[it]] <- dist1
      
      setTxtProgressBar(pb, it)
      
    }
    
    
    message("\n","Calculating PD ..." ,"\n")
    
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for(i in seq_len(iterations)){
      # compute PD
      pd_null <- suppressMessages(faith.pd(tree, iter.dist[[i]], method = pd$calculation.method))
      samples[, i] <- pd_null$pd.table$PD
      setTxtProgressBar(pb, i)
    }
  }
  
  
  cat("\n")

  ############################################################################################
  
  ##############################################################
  #                     NULL MODEL STATS                       #
  ##############################################################

  if (method == "global") {
    pd_obs <- sum(pd$pd.table$PD, na.rm = TRUE)
    null_pd <- as.vector(apply(samples, 2, sum, na.rm = TRUE))
    null_mean <- mean(null_pd)
    null_sd <- sd(null_pd)
  
    ## SAFE SES (avoid Inf when sd = 0)
  if (is.finite(null_sd) && null_sd == 0) {
    ses <- if (isTRUE(all.equal(pd_obs, null_mean))) 0 else NA_real_
  } else {
    ses <- (pd_obs - null_mean) / null_sd
  }
  
  ## Use only finite nulls in denominator
  n_eff <- sum(is.finite(null_pd))                                    # <<< ADDED
  if (n_eff == 0) {                                                   # <<< ADDED
    p_val <- NA_real_                                                 # <<< ADDED
  } else {                                                            # <<< ADDED
    p_val <- (sum(null_pd >= pd_obs, na.rm = TRUE) + 1) / (n_eff + 1) # <<< ADDED
  }                                                                   # <<< ADDED
  
  return(
    structure(
      list(pd.obs = pd_obs,
           null.pd = null_pd,
           SES = ses,
           Pvalue = p_val,
           Iterations = iterations,
           Iter.table = samples,
           Method = method),
      class = "vilma.null"
      )
    )
  }

  if (method == "cell") {
    pd_obs <- pd$pd.table$PD
    null_mean <- as.vector(apply(samples, 1, mean, na.rm = TRUE))
    null_sd <- as.vector(apply(samples, 1, sd,   na.rm = TRUE))
    ses <- (pd_obs - null_mean) / null_sd
  
    ## SAFE SES per-cell (avoid Inf when sd = 0)
    zero_sd <- is.finite(null_sd) & (null_sd == 0)
    same_as_mean <- abs(pd_obs - null_mean) < .Machine$double.eps^0.5
    ses[zero_sd &  same_as_mean] <- 0
    ses[zero_sd & !same_as_mean] <- NA_real_
  
    p_val <- numeric(nrow(samples))
  
  for (i in 1:nrow(samples)) {
    pd.obs1 <- pd_obs[i]
    null.pd <- samples[i, ]
    
    denom_i <- sum(is.finite(null.pd))                                # <<< ADDED
    if (denom_i == 0) {                                               # <<< ADDED
      p_val[i] <- NA_real_                                            # <<< ADDED
    } else {                                                          # <<< ADDED
      p_val[i] <- (sum(null.pd >= pd.obs1, na.rm = TRUE) + 1) / (denom_i + 1)  # <<< ADDED
    }                                                                 # <<< ADDED
  }
  
  cell.values <- data.frame(Cell = pd$pd.table$Cell,
                            PD = pd_obs,
                            NullPD = null_mean,
                            NullSD = null_sd,
                            SES = ses,
                            Pvalue = p_val)
  
  grid0 <- dist$grid
  values(grid0) <- NA
  set.values(grid0,
             as.numeric(as.character(cell.values$Cell)),
             as.numeric(ses))
  
  return(
    structure(
      list(CellValues = cell.values,
           Raster = grid0,
           Iterations = iterations,
           Iter.table = samples,
           Method = method),
      class = "vilma.null"
      )
    )
  }
  
}

