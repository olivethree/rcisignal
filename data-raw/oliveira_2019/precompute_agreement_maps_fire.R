## Precompute the 2x2 fire-palette agreement-maps figure for vignette
## section 12.x. Reads the cached per-trait signal matrices, draws one
## panel per trait via plot_agreement_map(palette = "fire"), composites
## on the base face, and writes a single 600 dpi PNG.
##
## Run from the rcisignal project root:
##   Rscript data-raw/oliveira_2019/precompute_agreement_maps_fire.R

suppressPackageStartupMessages(devtools::load_all("."))
library(png)

cache_dir <- file.path("vignettes", "precomputed", "oliveira_2019")
fig_dir   <- file.path("vignettes", "figures", "oliveira_2019")

stopifnot(dir.exists(cache_dir), dir.exists(fig_dir))

sm        <- readRDS(file.path(cache_dir, "signal_matrices.rds"))
base_face <- png::readPNG(file.path(fig_dir, "base_face.png"))
if (length(dim(base_face)) == 3L) {
  base_face <- 0.2126 * base_face[, , 1] +
               0.7152 * base_face[, , 2] +
               0.0722 * base_face[, , 3]
}

traits <- c("trust", "friendly", "competent", "dominant")
panel_titles <- c(
  trust     = "Trustworthy",
  friendly  = "Friendly",
  competent = "Competent",
  dominant  = "Dominant"
)

# Shared zlim across the four panels so the color scale is comparable.
t_maxes <- vapply(traits, function(tr) {
  out <- capture.output({
    res <- plot_agreement_map(sm[[tr]], palette = "fire",
                              main = "")
  })
  invisible(out)
  max(abs(res$t_map), na.rm = TRUE)
}, numeric(1L))
shared_zlim <- c(0, max(t_maxes))

out_path <- file.path(fig_dir, "agreement_maps_fire.png")
png(out_path, width = 8, height = 8, units = "in", res = 300, bg = "white")
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)

# Use layout() (not par(mfrow)) because plot_agreement_map() internally
# saves and restores par via par(no.readonly = TRUE); that restore would
# reset mfrow's panel counter on every call. layout() state survives
# par() restoration.
layout(matrix(1:4, nrow = 2, byrow = TRUE))

for (tr in traits) {
  plot_agreement_map(
    sm[[tr]],
    palette    = "fire",
    base_image = base_face,
    zlim       = shared_zlim,
    alpha_max  = 0.85,
    main       = panel_titles[[tr]]
  )
}

cat("Wrote ", out_path, "\n", sep = "")
