#### community_overlap_weighted() ####
# -------------------------------------------------------------------------
# A drop-in replacement for Ostats::community_overlap() that weights by an
# EXTERNAL abundance vector instead of the number of observations in each group.
#
# Why: the trait densities are estimated from an AUGMENTED sample (sparse species
# topped up to n = 20), so the observation count no longer reflects true
# abundance. community_overlap() uses that count in two places -- (a) the area
# under each species' density when normal = FALSE, and (b) the harmonic-mean
# weights on the pairwise overlaps -- so both are biased by augmentation. This
# version takes `abund`, a named vector (species -> true, effort-scaled abundance
# from NEON div data), and uses it in both places. Densities are still estimated
# from the (augmented) trait values; only the abundance weighting changes.
#
# It replicates community_overlap()'s machinery exactly (same grid limits, same
# stats::density defaults, same Sorensen overlap, same weightedMedian), verified
# to reproduce it to numerical precision when abund = the observation counts.
#
# ARGS
#   traits       numeric vector of trait values (e.g. log_dist_cm)
#   sp           species label per individual (same length as traits)
#   abund        named numeric vector: species -> abundance weight. Must cover
#                every eligible species (>= 2 individuals). Non-integer is fine.
#   normal       TRUE  -> each density area normalized to 1 (overlap_norm)
#                FALSE -> each density area proportional to abund[species] (overlap_unnorm)
#   output       "median" (default) or "mean"
#   weight_type  "hmean" (default), "mean", or "none" -- uses abund, not counts
#   randomize_weights  if TRUE, permute the pairwise weights (for null models)
#   density_args optional list; supports bw and n (defaults "nrd0", 512), matching Ostats
#
# RETURNS a single numeric O-statistic, or NA if < 2 eligible species.
# -------------------------------------------------------------------------

community_overlap_weighted <- function(traits, sp, abund,
                                       normal = TRUE, output = "median",
                                       weight_type = "hmean",
                                       randomize_weights = FALSE,
                                       density_args = list()) {

  traits <- as.numeric(traits)
  sp     <- as.character(sp)

  # clean: drop missing traits / labels, then species with < 2 individuals
  ok <- is.finite(traits) & !is.na(sp) & sp != ""
  traits <- traits[ok]; sp <- sp[ok]
  n_ind    <- table(sp)
  eligible <- names(n_ind)[n_ind > 1]
  keep     <- sp %in% eligible
  traits <- traits[keep]; sp <- sp[keep]

  uniquespp <- sort(unique(sp))
  nspp <- length(uniquespp)
  if (nspp < 2) return(NA)

  # abundance weight for each eligible species (must be supplied for all of them)
  w_sp <- abund[uniquespp]
  if (any(is.na(w_sp)))
    stop("community_overlap_weighted: no abundance for species ",
         paste(uniquespp[is.na(w_sp)], collapse = ", "))

  # common grid limits across all species (match community_overlap: extend range by +/- 0.5*range)
  rng  <- range(traits)
  grid <- rng + c(-0.5, 0.5) * diff(rng)

  bw <- if ("bw" %in% names(density_args)) density_args[["bw"]] else "nrd0"
  n  <- if ("n"  %in% names(density_args)) density_args[["n"]]  else 512

  # per-species density on the common grid; area = 1 (normal) or = true abundance (!normal)
  density_list <- lapply(uniquespp, function(s) {
    x <- traits[sp == s]
    d <- stats::density(x, from = grid[1], to = grid[2], bw = bw, n = n)
    y <- d$y
    if (!normal) y <- y * as.numeric(w_sp[s])   # area proportional to TRUE abundance
    list(x = d$x, y = y)
  })
  names(density_list) <- uniquespp

  # pairwise Sorensen overlaps + abundance weights
  combs    <- utils::combn(seq_len(nspp), 2)
  overlaps <- numeric(ncol(combs))
  wts      <- numeric(ncol(combs))
  for (k in seq_len(ncol(combs))) {
    a <- density_list[[combs[1, k]]]
    b <- density_list[[combs[2, k]]]
    wmin  <- pmin(a$y, b$y)
    total <- sfsmisc::integrate.xy(a$x, a$y) + sfsmisc::integrate.xy(b$x, b$y)
    inter <- sfsmisc::integrate.xy(a$x, wmin)
    overlaps[k] <- 2 * inter / total

    wa <- as.numeric(w_sp[uniquespp[combs[1, k]]])
    wb <- as.numeric(w_sp[uniquespp[combs[2, k]]])
    if (weight_type == "hmean")     wts[k] <- 2 / (1/wa + 1/wb)
    else if (weight_type == "mean") wts[k] <- wa + wb
    else                            wts[k] <- 1
  }

  if (randomize_weights) wts <- sample(wts)

  if (weight_type == "none")
    return(if (output == "median") stats::median(overlaps) else mean(overlaps))
  if (output == "median")
    return(matrixStats::weightedMedian(x = overlaps, w = wts))
  stats::weighted.mean(x = overlaps, w = wts)
}
