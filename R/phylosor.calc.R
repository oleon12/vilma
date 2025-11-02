#' Phylogenetic Sørensen Similarity (PhyloSor) between communities (cells)
#' @description
#' Computes pairwise PhyloSor similarity among grid cells and returns:
#' (i) a similarity matrix, (ii) its complementary dissimilarity matrix,
#' and (iii) rasters summarizing mean similarity/dissimilarity per cell
#' plus ordination (PCoA and NMDS) axes mapped to space.
#'
#' Two calculation modes are supported:
#' \itemize{
#'   \item \strong{Unweighted} (\code{abundance = FALSE}; default): presence–absence PhyloSor.
#'   \item \strong{Weighted} (\code{abundance = TRUE}): edge contributions are
#'         weighted by per-cell species \emph{abundances}. If \code{normalize = TRUE},
#'         abundances are converted to relative abundances within each cell
#'         (so each cell sums to 1), yielding a normalized, 0–1 scaled similarity.
#' }
#'
#' The PD convention for single-species cells follows \code{method}:
#' \code{"root"} uses root→tip distance; \code{"node"} uses the terminal branch;
#' \code{"exclude"} removes singletons from the analysis. For pairs sharing
#' exactly one species, behavior is controlled by \code{singleton_overlap}.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist} (e.g., from \code{points.to.raster()}),
#'   containing species occurrences per grid cell and a template raster.
#' @param method Character; PD convention for single-species cases:
#'   \itemize{
#'     \item \code{"root"} – singletons use root→tip; multi-species PD includes root→MRCA (default).
#'     \item \code{"node"} – singletons use terminal branch; multi-species PD is the minimal subtree.
#'     \item \code{"exclude"} – single-species cells are dropped before computing pairwise similarity.
#'   }
#' @param singleton_overlap Logical; when two cells share exactly one species:
#'   \itemize{
#'     \item \code{FALSE} – count that single shared species as shared PD under \code{method}.
#'     \item \code{TRUE} – treat that case as no shared PD (similarity = 0).
#'   }
#'   Note: \code{method="exclude"} is incompatible with \code{singleton_overlap = FALSE}.
#' @param abundance Logical; if \code{TRUE}, use abundance-weighted PhyloSor by
#'   partitioning shared/unique branch lengths with abundance weights (via
#'   \code{branch.partition.weighted()}). When \code{FALSE} (default), presence–absence is used.
#' @param normalize Logical; only relevant when \code{abundance = TRUE}. If \code{TRUE}
#'   (default), per-cell abundances are normalized to relative abundances before
#'   computing edge fractions. If \code{FALSE}, raw counts are used.
#'
#' @details
#' Let \eqn{PD_A}, \eqn{PD_B} be (convention-consistent) phylogenetic diversity
#' of communities A and B, and \eqn{PD_{shared}} their shared branch length.
#' The classic PhyloSor similarity is:
#' \deqn{\mathrm{PhyloSor}(A,B) = \frac{2\,PD_{shared}}{PD_A + PD_B}.}
#'
#' In the unweighted mode, presence–absence is used to compute \eqn{PD} and
#' \eqn{PD_{shared}}. In the weighted mode, for each tree edge \eqn{e} with length
#' \eqn{L_e}, the fractions of (relative) abundance descending from \eqn{e} in A and B
#' are computed (\eqn{p_A(e)}, \eqn{p_B(e)}). These yield three edge-weighted branch-length
#' sums: shared \eqn{a=\sum_e L_e \min[p_A(e),p_B(e)]}, unique-to-A
#' \eqn{b=\sum_e L_e \max[p_A(e)-p_B(e),0]}, and unique-to-B
#' \eqn{c=\sum_e L_e \max[p_B(e)-p_A(e),0]}. The weighted PhyloSor then uses
#' \eqn{PD_{shared}=a} and \eqn{PD_A+PD_B=2a+b+c}, so
#' \deqn{\mathrm{PhyloSor}_w(A,B)=\frac{2a}{2a+b+c}.}
#' With \code{normalize = TRUE}, \eqn{p_\cdot(e)} are computed from relative
#' abundances (summing to 1 within each cell), keeping the index in \eqn{[0,1]}.
#'
#' Only species present in both \code{tree} and \code{dist} are used; informative
#' messages are printed for mismatches. The similarity matrix has diagonal 1; the
#' dissimilarity matrix is \code{1 - similarity} with diagonal 0.
#'
#' After computing matrices, each cell is summarized by its mean similarity
#' (and \code{1 -} that value for mean dissimilarity). Ordinations (PCoA with
#' Cailliez correction, and NMDS) are done on the dissimilarity matrix; the
#' first two axes are rasterized to the study grid.
#'
#' @return An object of class \code{vilma.beta} with:
#' \itemize{
#'   \item \code{distribution} – Original species distribution (from \code{dist}).
#'   \item \code{similarity} – PhyloSor similarity as a \code{dist} object (labels = cell IDs).
#'   \item \code{dissimilarity} – \code{1 - similarity} as a \code{dist} object.
#'   \item \code{rasters} – A named list of \code{SpatRaster} layers:
#'         \code{mean.similarity}, \code{mean.dissimilarity},
#'         \code{pcoa.1}, \code{pcoa.2}, \code{ndms.1}, \code{ndms.2}.
#'   \item \code{pcoa.eig} – PCoA eigenvalues (first two axes).
#'   \item \code{ndms.stress} – NMDS stress.
#'   \item \code{calculation.method} – The resolved \code{method}.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item The tree must be rooted and have branch lengths.
#'   \item \code{abundance = TRUE} calls \code{branch.partition.weighted()}
#'         with \code{normalize} passed through; \code{abundance = FALSE}
#'         uses presence–absence via \code{branch.partition()}.
#'   \item \code{method = "exclude"} removes singletons prior to pairwise
#'         computations; such cells appear as \code{NA} in rasters.
#'   \item Ordinations are applied to dissimilarities; PCoA uses
#'         \code{ape::pcoa(..., correction = "cailliez")}.
#' }
#'
#' @references
#' Bryant JA, Lamanna C, Morlon H, Kerkhoff AJ, Enquist BJ, Green JL (2008)
#' Microbes on mountainsides: contrasting elevational patterns of bacterial and plant diversity.
#' \emph{Proceedings of the National Academy of Sciences} 105:11505–11511. (Introduces the PhyloSor concept/formulation.)
#'
#' Kembel SW, Cowan PD, Helmus MR, Cornwell WK, Morlon H, Ackerly DD, Blomberg SP, Webb CO (2010)
#' Picante: R tools for integrating phylogenies and ecology.
#' \emph{Bioinformatics} 26:1463–1464. (Implements \code{phylosor} and related functions.)
#'
#' Chiu C-H, Jost L, Chao A (2014)
#' Phylogenetic beta diversity, similarity, and differentiation measures based on Hill numbers.
#' \emph{Ecological Monographs} 84:21–44. (Abundance-based extensions and normalization considerations.)
#'
#' Gower JC (1966)
#' Some distance properties of latent root and vector methods used in multivariate analysis.
#' \emph{Biometrika} 53:325–338. (Classical basis of principal coordinates analysis.)
#'
#' Cailliez F (1983)
#' The analytical solution of the additive constant problem.
#' \emph{Psychometrika} 48:305–312. (PCoA negative-eigenvalue correction used here.)
#'
#' Kruskal JB (1964)
#' Nonmetric multidimensional scaling: a numerical method.
#' \emph{Psychometrika} 29:115–129. (Foundational NMDS procedure.)
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso
#' \code{\link{faith.pd}}, \code{\link{unifrac.calc}}, \code{\link{points_to_raster}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#'
#' # Presence–absence PhyloSor (default):
#' ps <- phylosor.calc(tree, dist)
#' print(ps)
#' plot(ps)
#' view.vilma(ps)
#'
#' # Abundance-weighted PhyloSor with relative abundances:
#' ps_w <- phylosor.calc(tree, dist, abundance = TRUE, normalize = TRUE)
#'
#' # Abundance-weighted without normalization (use raw counts):
#' ps_wc <- phylosor.calc(tree, dist, abundance = TRUE, normalize = FALSE)
#' }
#'
#'
#' @export

phylosor.calc <- function(tree, dist, method = c("root","node","exclude"), singleton_overlap = FALSE,
                          abundance = FALSE, normalize = TRUE){
  
  ##############################################################
  #                        VERIFICATION                        #
  ##############################################################
  
  if (!inherits(tree, "phylo")) {
    stop("Input 'tree' must be an object of class 'phylo'.")
  }
  
  if (is.null(tree$edge.length)) {
    stop("The phylogenetic tree must have branch lengths.")
  }
  
  if (!is.rooted(tree)) {
    stop("PhyloSor requires a rooted tree.")
  }
  
  if (!inherits(dist, "vilma.dist")) {
    stop("Input 'dist' must be an object of class 'vilma.dist'. See 'points_to_raster()' function.")
  }
  
  if(length(method)>2){
    method <- "root"
    message("Using method 'root' \n")
  }
  
  if(method == "exclude" && singleton_overlap == FALSE){
    stop("You can not use method exclude when singleton_overlap is FALSE. Only 'root' or 'node' are allowed.")
  }
  
  ##############################################################
  #                     DATA PREPARATION                       #
  ##############################################################
  
  ##### Dist Processing #####
  distM <- dist$distribution
  
  # 3. Find and use only species that match between tree and dist
  tree.species <- tree$tip.label
  dist.species <- unique(distM$Sp)
  
  missing.in.tree <- setdiff(dist.species, tree.species)
  missing.in.dist <- setdiff(tree.species, dist.species)
  
  if (length(missing.in.tree) > 0) {
    message("Note: ", length(missing.in.tree), 
            " species from 'dist' not found in 'tree': ", 
            paste(head(missing.in.tree, 5), collapse = ", "),
            ifelse(length(missing.in.tree) > 5, "...", ""))
  }
  
  if (length(missing.in.dist) > 0) {
    message("Note: ", length(missing.in.dist),
            " species from 'tree' not found in 'dist': ",
            paste(head(missing.in.dist, 5), collapse = ", "),
            ifelse(length(missing.in.dist) > 5, "...", ""))
  }
  
  
  ############################################################### 
  
  # Use only intersecting species
  common.species <- intersect(tree.species, dist.species)
  
  if (length(common.species) == 0) {
    stop("No species in common between the tree and distribution data.")
  }
  
  message("Using ", length(common.species), " species in common between tree and distribution.")
  distM <- distM[distM$Sp %in% common.species, ]
  
  ###############################################################
  
  # Handle cells with one species
  # two options: 1. remove them 2. root/node method
  
  if(method == "exclude"){
    
    countSp <- as.matrix(table(distM$Sp, distM$Cell))
    countSp <- apply(countSp, 2, function(x){length(x[x>0])})
    countSp <- names(countSp[which(countSp == 1)])
    
    distM <- distM[-which(distM$Cell %in% as.numeric(countSp)), ]
    dist0 <- dist
    dist0$distribution <- distM
    
  }else{
    dist0 <- dist
    dist0$distribution <- distM 
  }
  
  Cells <- sort(unique(dist0$distribution$Cell))
  if (length(Cells) < 2) stop("Not enough cells to compute PhyloSor after filtering.")
  
  #Create symetric matrix
  sym.matrix <- matrix(NA, nrow = length(unique(dist0$distribution$Cell)), 
                       ncol = length(unique(dist0$distribution$Cell)))
  
  rownames(sym.matrix) <- unique(dist0$distribution$Cell)
  colnames(sym.matrix) <- unique(dist0$distribution$Cell)
  
  #Get unique cells
  Cells <- unique(dist0$distribution$Cell)
  
  cat("\n")
  message("Calculating PhyloSor ...")
  cat("\n")
  
  desc_tips <- .desc_tips_by_edge(tree)
  
  pb <- txtProgressBar(min = 0, max = length(Cells) , style = 3)
  on.exit(close(pb), add = TRUE)
  
  for(a in seq_along(Cells)){
    
    #Get the distribution of community A
    DistA <- dist0
    DistA$distribution <- DistA$distribution[which(DistA$distribution$Cell%in%row.names(sym.matrix)[a]), ]
    
    #If community A with single species, then is the PD from tip to root/node (given single method)
    #If exclude method is selected, then also 
    #Here I use the whole original tree, if not, the distance to the root change
    if(length(unique(DistA$distribution$Sp)) == 1){
      PDA <- pd.taxon(tree = tree, sp = unique(DistA$distribution$Sp), method = method)
    }else{
      #More than one species, then calculate regular PD with community A's distribution
      PDA <- suppressMessages(faith.pd(tree = tree, dist = DistA, method = method))
      PDA <- PDA$pd.table$PD
    }
    
    for(b in seq_along(Cells)){
      
      if(a == b){
        sym.matrix[a,b] <- 1
        next
      }
      
      if(abundance == TRUE){
        
        # Build per-cell *counts* from distM (presence/absence -> counts = 1 each)
        commA <- distM$Sp[distM$Cell == Cells[a]]
        commB <- distM$Sp[distM$Cell == Cells[b]]
        
        # Named numeric abundance vectors (counts). If you ever pass true counts,
        # this code will automatically use them; with your data it’s all 1s.
        A_abund <- setNames(as.numeric(table(commA)), names(table(commA)))
        B_abund <- setNames(as.numeric(table(commB)), names(table(commB)))
        
        # Get a,b,c using your weighted partition (desc_tips precomputed)
        bl <- branch.partition.weighted(tree, desc_tips,
                                        commA_abund = A_abund,
                                        commB_abund = B_abund,
                                        normalize = normalize)
        
        a_len <- as.numeric(bl["blShare"])
        b_len <- as.numeric(bl["blA"])
        c_len <- as.numeric(bl["blB"])
        
        denom <- (2*a_len + b_len + c_len)
        sym.matrix[a,b] <- if (denom > 0) (2*a_len) / denom else 0
      }else{
        
        #Get distribution of community B
        DistB <- dist0
        DistB$distribution <- DistB$distribution[which(DistB$distribution$Cell%in%row.names(sym.matrix)[b]), ]
        
        #If community B with single species, then is the PD from tip to root/node (given single method)
        #Here I use the whole original tree, if not, the distance to the root change
        if(length(unique(DistB$distribution$Sp))==1){
          PDB <- pd.taxon(tree = tree, sp = unique(DistB$distribution$Sp), method = method)
        }else{
          #More than one species, then calculate regular PD with community B's distribution
          PDB <- suppressMessages(faith.pd(tree = tree, dist = DistB, method = method))
          PDB <- PDB$pd.table$PD
        }
        
        #Find common species
        common.sp <- intersect(unique(DistA$distribution$Sp), unique(DistB$distribution$Sp))
        
        if(length(common.sp) == 0){
          #If no common species, then phylosor is 0
          sym.matrix[a,b] <- 0
        }else{
          if(length(common.sp) == 1){
            
            if(singleton_overlap == FALSE){
              #If common species is 1, single cell species, then distance to tip to root/node (given single method)
              PDshare <- pd.taxon(tree = tree, sp = common.sp, method = method)
              #Calculate phylosor
              phylosor <- (2*PDshare) / (PDA + PDB)
              sym.matrix[a,b] <- phylosor 
            }else{
              sym.matrix[a,b] <- 0
            }
            
          }else{
            
            #If more than one common species, get distribution of only those species
            #Prune tree of common species
            TreeShare <- keep.tip(tree, common.sp)
            #Then, PDshare, will be the sum of branch lengths connecting those common species (Faith's PD)
            PDshare <- sum(TreeShare$edge.length)
            
            #Now, calculate phylosor
            phylosor <- (2*PDshare) / (PDA + PDB)
            sym.matrix[a,b] <- phylosor
          }
         }
       }
      }
    setTxtProgressBar(pb, a)
  }
  
  dis.matrix <- (1-sym.matrix)
  dis.matrix <- as.dist(dis.matrix)
  
  ################################################################
  #                        Mean PhyloSor                         #
  ################################################################
  
  mean.sym <- sym.matrix
  diag(mean.sym) <- NA
  
  mean.sym <- apply(mean.sym, 1, function(x){mean(x, na.rm =T)})
  mean.dis <- 1- mean.sym
  
  sym.matrix <- as.dist(sym.matrix)
  
  ################################################################
  #            Multidimensional reductions PCoA + NDMS           #
  ################################################################
  
  ############# Dissimilarity
  
  pcoa <- pcoa(dis.matrix, correction = "cailliez")
  
  eig.dis <- pcoa$values$Rel_corr_eig[1:2]
  
  pcoa.dis.1 <-pcoa$vectors[, 1]
  pcoa.dis.2 <-pcoa$vectors[, 2]
  
  ndms <- monoMDS(as.dist(dis.matrix), k = 2)
  ndms.stress.dis <- ndms$stress
  
  ndms.dis.1 <- ndms$points[, 1]
  ndms.dis.2 <- ndms$points[, 2]
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  #### Raster Mean PhyloSor####
  
  ### Sym
  mean.phylosor.sym <- dist$grid
  
  values(mean.phylosor.sym) <- NA
  
  set.values(mean.phylosor.sym, as.numeric(names(mean.sym)),  as.numeric(mean.sym))
  
  ### Dis
  mean.phylosor.dis <- dist$grid
  
  values(mean.phylosor.dis) <- NA
  
  set.values(mean.phylosor.dis, as.numeric(names(mean.dis)),  as.numeric(mean.dis))
  
  #### Raster PCoA ####
  
  pcoa.r1 <- dist$grid
  pcoa.r2 <- dist$grid
  
  values(pcoa.r1) <- NA
  values(pcoa.r2) <- NA
  
  set.values(pcoa.r1, as.numeric(names(pcoa.dis.1)),  as.numeric(pcoa.dis.1))
  set.values(pcoa.r2, as.numeric(names(pcoa.dis.2)),  as.numeric(pcoa.dis.2))
  
  #### Raster NDMS ####
  
  ndms.r1 <- dist$grid
  ndms.r2 <- dist$grid
  
  values(ndms.r1) <- NA
  values(ndms.r2) <- NA
  
  set.values(ndms.r1, as.numeric(names(ndms.dis.1)),  as.numeric(ndms.dis.1))
  set.values(ndms.r2, as.numeric(names(ndms.dis.2)),  as.numeric(ndms.dis.2))
  
  ##############################################################
  #                          Output                            #
  ##############################################################
  
  structure(
    list(
      distribution = dist$distribution,
      similarity = sym.matrix,
      dissimilarity = dis.matrix,
      rasters = list(mean.similarity = mean.phylosor.sym,
                     mean.dissimilarity = mean.phylosor.dis,
                     pcoa.1 = pcoa.r1,
                     pcoa.2 = pcoa.r2,
                     ndms.1 = ndms.r1,
                     ndms.2 = ndms.r2),
      pcoa.eig = eig.dis,
      ndms.stress =ndms.stress.dis,
      calculation.method = method,
      algorithm = "PhyloSor"
    ),
    class = "vilma.beta"
  )
}
