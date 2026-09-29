#### PHYLOGENETIC DISTANCE PER COMMUNITY ####
# -------------------------------------------------------------------------
# Phylogenetic metrics for every plot_year and site_year community.
# PD and its SES come from picante::ses.pd; MPD and MNTD are observed only. Community definition matches script 6: species observed in
# the focal unit that year with a true (effort-scaled) abundance.
#
# METRICS (presence-based):
#   pd      Faith's PD without the root path (picante::ses.pd, pd.obs)
#   pd_ses  standardized effect size of PD vs the null (pd.obs.z)
#   pd_p    quantile of observed PD in the null (pd.obs.p)
#   mpd   mean pairwise phylogenetic distance
#   mntd  mean nearest-taxon phylogenetic distance
# All three are tip-to-tip path lengths, so they do not depend on rooting.
# Branch lengths are substitutions (tree is not ultrametric).
#
# OUTPUTS:
#   ./Outputs/plot_PhyloDistance.csv
#   ./Outputs/site_PhyloDistance.csv
# -------------------------------------------------------------------------

library(ape)
library(picante)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

#### 0. SETTINGS ####
TREE_FILE <- "./Data/Tree.DNA.droppedseqs.newick"
YEARS     <- c(2018, 2019)
LEVELS    <- c("plot", "site")

NULL_MODEL <- "taxa.labels"   # shuffle tips across the (pruned) tree
RUNS       <- 999
SEED       <- 517


#### 1. READ TREE AND HARMONIZE TIP LABELS ####
# Tips are "Genus_species_N" (N = trailing sequence index). Strip the number,
# then underscores -> spaces.
tree <- read.tree(TREE_FILE)
tree$tip.label <- gsub("_", " ", sub("_\\d+$", "", tree$tip.label))

# drop unnamed tips (one tip is labelled only "_138")
bad_tip <- !grepl("^[A-Z][a-z]+ [a-z]", tree$tip.label)
if (any(bad_tip)) tree <- drop.tip(tree, which(bad_tip))

# one tip per species
dup <- duplicated(tree$tip.label)
if (any(dup)) tree <- drop.tip(tree, which(dup))

cophen <- cophenetic(tree)


#### 2. READ CLEAN SPECIMEN DATA ####
all_elytra <- read.csv("./Data/BodysizeCombinedClean.csv")
all_elytra <- subset(all_elytra, yearCollected %in% YEARS)
all_elytra <- all_elytra[!is.na(all_elytra$scientificName_Species) &
                           all_elytra$scientificName_Species != "" &
                           is.finite(all_elytra$cm_elytra_max_length), ]

name_fix <- c("Amara crassipina"   = "Amara crassispina")   # only if not Puerto Rico
hit <- all_elytra$scientificName_Species %in% names(name_fix)
all_elytra$scientificName_Species[hit] <- name_fix[all_elytra$scientificName_Species[hit]]

# species in the data but not in the tree (check for synonymy / spelling)
missing_sp <- setdiff(sort(unique(all_elytra$scientificName_Species)), tree$tip.label)
message(length(missing_sp), " species missing from tree")
missing_sp
write.csv(data.frame(scientificName_Species = missing_sp),
          "./Outputs/PhyloMissingSpecies.csv", row.names = FALSE)


#### 3. METRIC HELPER ####
phylo_metrics <- function(sp) {
  out <- c(mpd = NA_real_, mntd = NA_real_)
  if (length(sp) < 2) return(out)
  d <- cophen[sp, sp]
  out["mpd"]  <- mean(d[upper.tri(d)])
  diag(d)     <- Inf
  out["mntd"] <- mean(apply(d, 1, min))
  out
}


#### 4. LOOP OVER LEVEL x YEAR, ONE CSV PER LEVEL ####
for (LEVEL in LEVELS) {
  FOCAL_COL <- if (LEVEL == "plot") "plotID" else "siteID"
  rows <- list()
  comm_long <- list()
  
  for (yr in YEARS) {
    dat <- unique(all_elytra[all_elytra$yearCollected == yr,
                             c(FOCAL_COL, "scientificName_Species")])
    
    # community = species with a true abundance in that unit-year
    abund_tab <- read.csv(sprintf("./Data/%s_abund_%d.csv", LEVEL, yr))
    abund_key <- paste(abund_tab[[FOCAL_COL]], abund_tab$scientificName_Species, sep = "|")
    dat <- dat[paste(dat[[FOCAL_COL]], dat$scientificName_Species, sep = "|") %in% abund_key, ]
    
    for (f in sort(unique(dat[[FOCAL_COL]]))) {
      sp_all  <- unique(dat$scientificName_Species[dat[[FOCAL_COL]] == f])
      sp_tree <- intersect(sp_all, tree$tip.label)
      m <- phylo_metrics(sp_tree)
      if (length(sp_tree)) comm_long[[length(comm_long) + 1]] <-
        data.frame(Assemblage = paste0(f, "_", yr), sp = sp_tree)
      rows[[length(rows) + 1]] <- data.frame(
        focal = f, Year = yr, Assemblage = paste0(f, "_", yr),
        n_sp = length(sp_all), n_sp_tree = length(sp_tree),
        mpd = m["mpd"], mntd = m["mntd"],
        row.names = NULL)
    }
  }
  
  out <- do.call(rbind, rows)
  
  # ---- PD and SES of PD via picante::ses.pd ----
  # Presence/absence matrix: rows = Assemblage, cols = species in the tree.
  # The tree is pruned to species present in this matrix, so taxa.labels
  # shuffles across the dataset's species rather than all 654 tips.
  # include.root = FALSE because the tree is unrooted.
  cl   <- do.call(rbind, comm_long)
  comm <- as.matrix(unclass(table(cl$Assemblage, cl$sp)))
  comm[comm > 0] <- 1
  tree_lvl <- keep.tip(tree, colnames(comm))
  comm <- comm[, tree_lvl$tip.label]
  
  set.seed(SEED)
  sp_res <- ses.pd(comm, tree_lvl, null.model = NULL_MODEL,
                   runs = RUNS, include.root = FALSE)
  sp_res <- data.frame(Assemblage = rownames(sp_res),
                       pd = sp_res$pd.obs, pd_ses = sp_res$pd.obs.z,
                       pd_p = sp_res$pd.obs.p)
  out <- merge(out, sp_res, by = "Assemblage", all.x = TRUE)
  
  # PD and its SES are undefined for fewer than two species in the tree
  out[out$n_sp_tree < 2, c("pd", "pd_ses", "pd_p")] <- NA
  names(out)[names(out) == "focal"] <- FOCAL_COL
  write.csv(out, sprintf("./Outputs/%s_PhyloDistance.csv", LEVEL), row.names = FALSE)
  message("wrote ", LEVEL, "_PhyloDistance.csv (", nrow(out), " rows)")
}
