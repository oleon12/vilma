# inst/tinytest/test_vilma.R

if (!requireNamespace("tinytest", quietly = TRUE)) {
  stop("Package 'tinytest' is required to run these tests.")
}

library(tinytest)
library(ape)
library(terra)
library(sf)

#############
# Fake data #
#############

# Fake phylogeny
fake_tree <- read.tree(text = "((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);")

fake_occ <- data.frame( Sp = c( "sp1", "sp1", "sp1","sp2", "sp2", "sp2","sp3", "sp3", "sp3","sp4", "sp4", "sp4"),
                        Lon = c(0.2, 0.8, 1.4, 1.2, 1.8, 0.4, 0.3, 1.3, 1.7, 0.5, 1.5, 0.9),
                        Lat = c(  0.2, 0.8, 1.6, 1.2, 1.8, 0.4,1.5, 0.7, 1.3,0.5, 1.5, 1.1) )

# Fake polygons
p1 <- st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))
p2 <- st_polygon(list(rbind(c(1, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 0))))
p3 <- st_polygon(list(rbind(c(0, 1), c(1, 1), c(1, 2), c(0, 2), c(0, 1))))
p4 <- st_polygon(list(rbind(c(1, 1), c(2, 1), c(2, 2), c(1, 2), c(1, 1))))

fake_pols <- st_sf(
  species = c("sp1", "sp2", "sp3", "sp4"),
  geometry = st_sfc(p1, p2, p3, p4),
  crs = 4326
)

#####################################
# Tests for fake objects themselves #
#####################################

expect_true(inherits(fake_tree, "phylo"))
expect_equal(length(fake_tree$tip.label), 4)
expect_equal(sort(fake_tree$tip.label), c("sp1", "sp2", "sp3", "sp4"))

expect_true(is.data.frame(fake_occ))
expect_equal(colnames(fake_occ), c("Sp", "Lon", "Lat"))
expect_equal(nrow(fake_occ), 12)

expect_true(all(unique(fake_occ$Sp) %in% fake_tree$tip.label))

expect_true(inherits(fake_pols, "sf"))
expect_equal(nrow(fake_pols), 4)


####################
# Points to raster #
####################

v_dist <- points_to_raster(points = fake_occ, res = 0.5)
expect_true(inherits(v_dist, "vilma.dist"))

#################
# Alpha indices #
#################

#Faith PD
faith_pd <- faith.pd(tree = fake_tree, dist = v_dist, method = "root")
expect_true(inherits(faith_pd, "vilma.pd"))

#MPD
mpd <- mpd.calc(tree = fake_tree, dist = v_dist, method = "root", abundance = F)
expect_true(inherits(mpd, "vilma.pd"))

#MNTD
mntd <- mntd.calc(tree = fake_tree, dist = v_dist, method = "root", abundance = F)
expect_true(inherits(mntd, "vilma.pd"))

#RaoQ
rao <- rao.calc(tree = fake_tree, dist = v_dist, abundance = F, scale01 = TRUE)
expect_true(inherits(rao, "vilma.pd"))

#Phylogenetic Endemism
pe <- pe.calc(tree = fake_tree, dist = v_dist, RPE = T, faith.method = "root")
expect_true(inherits(pe, "vilma.pd"))

#NTI
nti <- nti.calc(tree = fake_tree, dist = v_dist, iterations = 5, abundance = F)
expect_true(inherits(nti, "vilma.pd"))

#NRI
nri <- nri.calc(tree = fake_tree, dist = v_dist, iterations = 5, abundance = F)
expect_true(inherits(nri, "vilma.pd"))

###############
# Null Models #
###############

#Faith PD
faith_null <- faith.pd.null(pd = faith_pd, tree = fake_tree, dist = v_dist, iterations = 5, method = "global", sampling = "taxa.label")
expect_true(inherits(faith_null, "vilma.null"))

#MPD
mpd_null <- mpd.calc.null(pd = mpd, tree = fake_tree, dist = v_dist, iterations = 5, method = "global", sampling = "taxa.label")
expect_true(inherits(mpd_null, "vilma.null"))

#MNTD
mntd_null <- mntd.calc.null(pd = mntd, tree = fake_tree, dist = v_dist, iterations = 5, method = "global", sampling = "taxa.label")
expect_true(inherits(mntd_null, "vilma.null"))

#RaoQ
rao_null <- rao.calc.null(pd = rao, tree = fake_tree, dist = v_dist, iterations = 5, method = "global", sampling = "taxa.label", abundance = F, scale01 = T)
expect_true(inherits(rao_null, "vilma.null"))

#Phylogenetic Endemism
pe_null <- pe.calc.null(pd = pe, tree = fake_tree, dist = v_dist, iterations = 5, method = "global", sampling = "taxa.label")
expect_true(inherits(pe_null, "vilma.null"))

################
# Beta indices #
################

#Phylosor
phylosor <- phylosor.calc(tree = fake_tree, dist = v_dist)
expect_true(inherits(phylosor, "vilma.beta"))

#Phylobeta
phyb <- phylo.beta(tree = fake_tree, dist = v_dist)
expect_true(inherits(phyb, "vilma.beta"))

#Beta-MPD
mpdb <- beta.mpd(tree = fake_tree, dist = v_dist)
expect_true(inherits(mpdb, "vilma.beta"))

#Beta-MNTD
mntdb <- beta.mntd(tree = fake_tree, dist = v_dist)
expect_true(inherits(mntdb,"vilma.beta"))

#Beta-RaoQ
raob <- rao.beta(tree = fake_tree, dist = v_dist)
expect_true(inherits(raob,"vilma.beta"))

#UniFrac
unib <- unifrac.calc(tree = fake_tree, dist = v_dist)
expect_true(inherits(unib,"vilma.beta"))

###################
# Regionalization #
###################

region <- vilma.regionalize(beta = phylosor, k = 2, method = "network", optimize = F)
expect_true(inherits(region, "vilma.region"))
