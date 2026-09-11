#### community_depth() ####
# -------------------------------------------------------------------------
# Community-wide "overlap depth": at each point on the trait axis, how many
# species co-occupy it, averaged over the occupied range. A community-wide
# alternative to the pairwise O-statistic that does NOT floor out or become
# zero-inflated at high richness -- the pairwise mean divides the overlap
# signal by S(S-1)/2 pairs, whereas this divides one summed curve by the
# occupied width.
#
# Construction (peak-scaled depth):
#   1. estimate one density per species on a COMMON grid -- identical grid /
#      bw / n to community_overlap_weighted() so the two metrics see the same
#      densities;
#   2. peak-scale each density to its own maximum, so every species contributes
#      at most 1 where it is densest and tapers to 0 in its tails
#      (o_i(x) = f_i(x) / max f_i);
#   3. sum across species  ->  D(x) = number of species co-occupying x;
#   4. average D(x) over the occupied trait range (2.5-97.5% span of the pooled
#      individuals -- the SAME support as niche_range).
#
# Each species contributes equally (presence-based): depth is deliberately NOT
# abundance-weighted, which keeps it robust to the pitfall activity-density bias
# that contaminates the abundance weights used by the overlap metrics. For an
# abundance-weighted variant, scale o_i by w_sp[s] / max(w_sp) before summing.
#
# INTERPRETATION: observed depth scales with richness by construction (more
# species in a bounded axis -> more pile-up), so the quantity to use downstream
# is the SES against the POOL null, which holds richness constant. Depth HIGHER
# than the same-richness pool = trait space more packed than random
# (clustering / filtering-consistent); LOWER = more even (overdispersion /
# limiting-similarity-consistent). The raw observed value is descriptive only.
#
# ARGS
#   traits        numeric vector of trait values (e.g. log_dist_cm)
#   sp            species label per individual (same length as traits)
#   occ_probs     quantiles defining the occupied range (default 2.5-97.5%,
#                 matching niche_range)
#   density_args  optional list; supports bw and n (defaults "nrd0", 512),
#                 matching community_overlap_weighted / Ostats
#
# RETURNS a single numeric depth (mean species co-occupancy over the occupied
# range), or NA if < 2 eligible species.
# -------------------------------------------------------------------------

community_depth <- function(traits, sp,
                            occ_probs = c(0.025, 0.975),
                            density_args = list()) {

  traits <- as.numeric(traits)
  sp     <- as.character(sp)

  # clean: drop missing traits / labels, then species with < 2 individuals
  # (identical eligibility rule to community_overlap_weighted)
  ok <- is.finite(traits) & !is.na(sp) & sp != ""
  traits <- traits[ok]; sp <- sp[ok]
  n_ind    <- table(sp)
  eligible <- names(n_ind)[n_ind > 1]
  keep     <- sp %in% eligible
  traits <- traits[keep]; sp <- sp[keep]

  uniquespp <- sort(unique(sp))
  if (length(uniquespp) < 2) return(NA)

  # common grid limits across all species (identical to community_overlap_weighted:
  # extend the data range by +/- 0.5*range). from/to/n are fixed, so every
  # species' density is evaluated on the SAME x grid and the y's are summable.
  rng  <- range(traits)
  grid <- rng + c(-0.5, 0.5) * diff(rng)
  bw <- if ("bw" %in% names(density_args)) density_args[["bw"]] else "nrd0"
  n  <- if ("n"  %in% names(density_args)) density_args[["n"]]  else 512

  # peak-scaled occupancy summed across species -> D(x)
  xg    <- NULL
  depth <- numeric(n)
  for (s in uniquespp) {
    d  <- stats::density(traits[sp == s], from = grid[1], to = grid[2], bw = bw, n = n)
    if (is.null(xg)) xg <- d$x
    mx <- max(d$y)
    if (mx > 0) depth <- depth + d$y / mx          # o_i(x) in [0, 1]
  }

  # average D(x) over the occupied trait range (same support as niche_range)
  occ      <- quantile(traits, occ_probs)
  in_range <- xg >= occ[1] & xg <= occ[2]
  if (!any(in_range)) return(NA)
  mean(depth[in_range])
}
