#' Correlate phylogenetic diversity with species richness in a \code{vilma.pd} object
#'
#' Fits simple statistical models relating phylogenetic diversity (PD) to species
#' richness (SR) using the per-cell summary table stored in a \code{vilma.pd}
#' object. The function fits three model families—linear model (LM), generalized
#' linear model (GLM), and generalized additive model (GAM)—and selects the best
#' model within the GLM family and within the GAM family using AIC. It then
#' compares the best LM, best GLM, and best GAM using AIC and returns a
#' \code{vilma.cor} object for downstream printing/plotting.
#'
#' @param vilma.pd A \code{vilma.pd} object containing \code{vilma.pd$pd.table} with
#'   at least two numeric columns: \code{PD} (phylogenetic diversity) and \code{SR}
#'   (species richness).
#'
#' @details
#' The following models are fit:
#' \itemize{
#'   \item \strong{LM:} \code{lm(PD ~ SR)}.
#'   \item \strong{GLM:} Gamma(log), inverse Gaussian(log), and Gaussian(log); the
#'   GLM with minimum AIC is retained.
#'   \item \strong{GAM:} a linear term \code{gam(PD ~ SR)} and a smooth term
#'   \code{gam(PD ~ s(SR))}; the GAM with minimum AIC is retained.
#' }
#' The final AIC comparison among \code{lm}, selected \code{glm}, and selected
#' \code{gam} is returned in \code{$AIC}.
#'
#' @return An object of class \code{vilma.cor}, a list with:
#' \describe{
#'   \item{Data}{A data frame with columns \code{pd.vals} and \code{sr.vals} used for fitting.}
#'   \item{Models}{A list with elements \code{lm}, \code{glm}, and \code{gam} containing
#'     the fitted model objects.}
#'   \item{AIC}{Named numeric vector with AIC values for \code{lm}, \code{glm}, and \code{gam}.}
#' }
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#'   \code{\link[stats]{AIC}}, \code{\link[mgcv]{gam}},
#'   \code{\link{plot.vilma.cor}}
#'
#' @examples
#' \dontrun{
#' cor_out <- vilma.cor(vilma_pd)
#' cor_out$AIC
#' plot(cor_out)
#' }
#'
#' @export

vilma.cor <- function(vilma.pd){
  if(!inherits(vilma.pd, "vilma.pd")){
    stop("Input is not a 'vilma.pd' object")
  }

  pd.vals <- vilma.pd$pd.table$PD
  sr.vals <- vilma.pd$pd.table$SR

  lm_out <- lm(pd.vals ~ sr.vals)
  glm_outLOG <- glm(pd.vals ~ sr.vals, family = Gamma(link = "log"))
  glm_outIG <- glm(pd.vals ~ sr.vals, family = inverse.gaussian(link = "log"))
  glm_outGau <- glm(pd.vals ~ sr.vals, family = gaussian(link = "log"))

  glm_tmp <- list(glm_outLOG, glm_outIG, glm_outGau)
  glm_aic <- c(LOG = AIC(glm_outLOG), INV = AIC(glm_outIG), GAU = AIC(glm_outGau))
  glm_out <- glm_tmp[[which(glm_aic == min(glm_aic))]]

  gam_n <- gam(pd.vals ~ sr.vals)
  gam_s <- gam(pd.vals ~ s(sr.vals))

  gam_tmp <- list(gam_n, gam_s)
  gam_aic <- c(Normal = AIC(gam_n), Smooth = AIC(gam_s))
  gam_out <- gam_tmp[[which(gam_aic == min(gam_aic))]]

  out <- list(lm = lm_out,
              glm = glm_out,
              gam = gam_out)

  AIC_vals <- c(lm = AIC(lm_out), glm = AIC(glm_out), gam = AIC(gam_out))

  stats_out <- list(Data = data.frame(pd.vals = pd.vals, sr.vals = sr.vals),
                    Models = out,
                    AIC = AIC_vals)

  return(structure(stats_out, class = "vilma.cor"))
}
