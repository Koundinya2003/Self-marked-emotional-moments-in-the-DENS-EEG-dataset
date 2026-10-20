# Region definitions on the EGI HydroCel GSN 128 net (template coordinates
# stored in the .set files; cm, +X towards nasion, +Y towards the left ear).
#
# Each region is the union of anchor channels (approximate 10-20 equivalents
# for this net) and every scalp channel within REGION_RADIUS_CM of an anchor.
# The outer ring of face, cheek and neck channels is never used for features
# because it is dominated by muscle and ocular activity.

OUTER_RING <- paste0("E", c(1, 8, 14, 17, 21, 25, 32, 38, 43, 44, 48, 49, 56, 57, 63, 64, 68, 69,
                            73, 74, 81, 82, 88, 89, 94, 95, 99, 100, 107, 113, 114, 119, 120,
                            121, 125, 126, 127, 128))

# Periocular channels used by the ICA eye-component rule and for EOG proxies
EOG_CHANS <- list(
  veog_left = c("E25", "E127"),   # above and below the left eye
  veog_right = c("E8", "E126"),   # above and below the right eye
  heog = c("E128", "E125")        # left and right outer canthi
)

REGION_ANCHORS <- list(
  frontal   = c(F3 = "E24", Fz = "E11", F4 = "E124"),
  central   = c(C3 = "E36", Cz = "E129", C4 = "E104"),
  posterior = c(P3 = "E52", Pz = "E62", P4 = "E92", O1 = "E70", Oz = "E75", O2 = "E83"),
  left      = c(F7 = "E33", T7 = "E45", P7 = "E58"),
  right     = c(F8 = "E122", T8 = "E108", P8 = "E96"),
  # frontal alpha asymmetry clusters
  faa_left  = c(F3 = "E24"),
  faa_right = c(F4 = "E124")
)
REGION_RADIUS_CM <- 3.0

define_regions <- function(chans, radius = REGION_RADIUS_CM) {
  scalp <- chans[chans$type == "EEG" & !(chans$label %in% OUTER_RING), ]
  xyz <- as.matrix(scalp[, c("X", "Y", "Z")])
  rownames(xyz) <- scalp$label
  lapply(REGION_ANCHORS, function(a) {
    d <- sapply(a, function(ch) sqrt(colSums((t(xyz) - xyz[ch, ])^2)))
    sort_ch(rownames(xyz)[apply(d, 1, min) <= radius])
  })
}

sort_ch <- function(x) x[order(as.integer(sub("E", "", x)))]

# Channel lists written into the preregistration (section 4.4). The pipeline
# stops if the geometric rule ever yields anything else.
PREREG_REGIONS <- list(
  frontal   = paste0("E", c(3, 4, 11, 16, 19, 20, 23, 24, 27, 28, 117, 118, 123, 124)),
  central   = paste0("E", c(7, 29, 30, 31, 35, 36, 41, 55, 80, 103, 104, 105, 106, 110, 111, 129)),
  posterior = paste0("E", c(42, 47, 51, 52, 60, 62, 65, 67, 70, 72, 75, 77, 83, 85, 90, 92, 93, 97, 98)),
  left      = paste0("E", c(33, 45, 50, 58, 65)),
  right     = paste0("E", c(90, 96, 101, 108, 122)),
  faa_left  = paste0("E", c(19, 20, 23, 24, 27, 28)),
  faa_right = paste0("E", c(3, 4, 117, 118, 123, 124))
)
check_regions <- function(regions) {
  stopifnot(identical(names(regions), names(PREREG_REGIONS)),
            all(mapply(identical, regions, PREREG_REGIONS)))
  invisible(TRUE)
}
