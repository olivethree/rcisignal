## Precompute the trait-CI correlation matrix and the ICC(3,k)
## vs group-mean infoVal-z scatter shown in vignette sections 7.2
## and 12.6.1. Reads cached signal_matrices.rds, infoval_full.rds,
## reliability_table.rds; writes:
##
##   vignettes/precomputed/oliveira_2019/trait_ci_correlations.rds
##   vignettes/precomputed/oliveira_2019/icc_vs_groupz.rds
##   vignettes/figures/oliveira_2019/trait_ci_correlogram.png
##   vignettes/figures/oliveira_2019/icc_vs_groupz.png

suppressPackageStartupMessages({
  library(devtools)
})
devtools::load_all(".")

fig_dir   <- "vignettes/figures/oliveira_2019"
cache_dir <- "vignettes/precomputed/oliveira_2019"

dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Trait CI correlation matrix --------------------------------

message("Computing pairwise correlations across the ten group-mean CIs ...")

sm <- readRDS(file.path(cache_dir, "signal_matrices.rds"))
group_cis <- sapply(sm, rowMeans)               # n_pixels x n_traits

# Match Oliveira et al. (2019, p. 892): CIs were masked with an oval
# over face / hair / ears before being correlated. Off-face pixels
# carry sinusoid energy from the noise pool and dilute every pair if
# included; restricting to the oval recovers the paper's Table 2
# values to +/- 0.02 across 17 of 18 spot-checked cells.
oval      <- make_face_mask(c(256L, 256L), region = "full")
trait_cor <- cor(group_cis[oval, , drop = FALSE])

saveRDS(trait_cor, file.path(cache_dir, "trait_ci_correlations.rds"))

# Restrict the displayed correlogram to the four traits used in
# the vignette's worked-example contrasts (trust, friendly,
# dominant, competent). Order chosen so positive blocks cluster
# top-left.
trait_order <- c("trust", "friendly", "competent", "dominant")
stopifnot(all(trait_order %in% colnames(trait_cor)))
trait_cor <- trait_cor[trait_order, trait_order]

# ---- 2. Upper-triangle correlogram PNG -----------------------------

message("Rendering trait-CI correlogram ...")

correlogram_upper <- function(M, out_path,
                              palette = grDevices::hcl.colors(
                                256L, "RdBu")) {
  n <- nrow(M)
  display <- M
  display[lower.tri(display, diag = TRUE)] <- NA_real_

  grDevices::png(out_path, width = 1200, height = 1200, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(0.5, 7.5, 7.5, 6.5) + 0.1)

  graphics::image(
    1:n, 1:n, t(display[n:1, ]),
    col       = palette, zlim = c(-1, 1),
    axes      = FALSE, xlab = "", ylab = "",
    asp       = 1, useRaster = TRUE
  )

  # Numbers in upper-triangle cells.
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      v <- display[i, j]
      if (!is.na(v)) {
        graphics::text(
          x = j, y = n - i + 1L,
          labels = sprintf("%+.2f", v),
          col = if (abs(v) > 0.6) "white" else "grey20",
          cex = 0.8
        )
      }
    }
  }

  # Trait labels.
  graphics::axis(3, at = seq_len(n), labels = colnames(M),
                 las = 2, tick = FALSE, line = -0.5, cex.axis = 0.85)
  graphics::axis(2, at = rev(seq_len(n)), labels = rownames(M),
                 las = 1, tick = FALSE, line = -0.5, cex.axis = 0.85)

  # Colour bar in right margin.
  usr <- graphics::par("usr")
  x_left  <- usr[2] + (usr[2] - usr[1]) * 0.02
  x_right <- usr[2] + (usr[2] - usr[1]) * 0.06
  y_bot   <- usr[3]; y_top <- usr[4]
  ys      <- seq(y_bot, y_top, length.out = length(palette) + 1L)
  graphics::par(xpd = TRUE)
  for (i in seq_along(palette)) {
    graphics::rect(x_left, ys[i], x_right, ys[i + 1L],
                   col = palette[i], border = NA)
  }
  graphics::rect(x_left, y_bot, x_right, y_top,
                 border = "grey60", lwd = 0.5)
  ticks <- seq(-1, 1, by = 0.5)
  tick_y <- y_bot + (ticks + 1) / 2 * (y_top - y_bot)
  graphics::segments(x_right, tick_y,
                     x_right + (usr[2] - usr[1]) * 0.01,
                     tick_y, col = "grey60", lwd = 0.5)
  graphics::text(x_right + (usr[2] - usr[1]) * 0.015, tick_y,
                 labels = format(ticks), pos = 4, cex = 0.7,
                 col = "grey20")
  graphics::par(xpd = FALSE)
}

correlogram_upper(trait_cor,
                  file.path(fig_dir, "trait_ci_correlogram.png"))

# ---- 3. ICC(3,k) vs group-mean z lm fit ----------------------------

message("Fitting linear model: group_z ~ ICC(3,k) across the 10 traits ...")

iv  <- readRDS(file.path(cache_dir, "infoval_full.rds"))
rel <- readRDS(file.path(cache_dir, "reliability_table.rds"))
df  <- merge(rel, iv, by = "trait")
df  <- df[, c("trait", "icc_3_k", "group_z", "median_z")]

ct  <- stats::cor.test(df$icc_3_k, df$group_z)
fit <- stats::lm(group_z ~ icc_3_k, data = df)

icc_grid <- seq(min(df$icc_3_k) - 0.05,
                max(df$icc_3_k) + 0.05, length.out = 200L)
pred <- as.data.frame(stats::predict(
  fit, newdata = data.frame(icc_3_k = icc_grid),
  interval = "confidence", level = 0.95
))
pred$icc_3_k <- icc_grid

icc_vs_groupz <- list(
  data    = df,
  fit     = fit,
  pred    = pred,
  cor     = c(estimate  = unname(ct$estimate),
              ci_lo     = ct$conf.int[1],
              ci_hi     = ct$conf.int[2],
              t         = unname(ct$statistic),
              df        = unname(ct$parameter),
              p_value   = ct$p.value)
)
saveRDS(icc_vs_groupz, file.path(cache_dir, "icc_vs_groupz.rds"))

# ---- 4. Scatter PNG with fit + 95% CI band -------------------------

message("Rendering ICC(3,k) vs group-mean z scatter ...")

scatter_icc_vs_groupz <- function(d, pred, cor_stats, out_path) {
  grDevices::png(out_path, width = 1200, height = 950, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(4.4, 4.6, 1.2, 1.2) + 0.1)

  # Extra room on the left so the bottom-cluster labels don't clip
  # the y-axis; the fit line/band still only spans the data range.
  xlim <- c(min(d$icc_3_k) - 0.12, max(d$icc_3_k) + 0.07)
  ylim <- range(c(pred$lwr, pred$upr, d$group_z))
  ylim[1] <- min(ylim[1], 0)

  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "ICC(3,k) at the group level",
       ylab = "Group-mean infoVal z",
       las = 1, bty = "n",
       panel.first = graphics::grid(nx = NULL, ny = NULL,
                                    col = "grey92", lty = 1, lwd = 0.5))

  # 95 % CI band.
  graphics::polygon(
    x = c(pred$icc_3_k, rev(pred$icc_3_k)),
    y = c(pred$lwr,     rev(pred$upr)),
    col = grDevices::adjustcolor("steelblue", alpha.f = 0.20),
    border = NA
  )

  # Fit line.
  graphics::lines(pred$icc_3_k, pred$fit, col = "steelblue3", lwd = 2)

  # Reference lines: group_z = 1.96 (1-tail z threshold) and 0.
  graphics::abline(h = 1.96, col = "grey70", lty = 3, lwd = 0.7)
  graphics::abline(h = 0,    col = "grey80", lty = 1, lwd = 0.5)

  # Points.
  graphics::points(d$icc_3_k, d$group_z, pch = 21,
                   bg = "white", col = "grey25", lwd = 1.2, cex = 1.4)

  # Trait labels. Hand-placed so the bottom-left cluster stays
  # readable (incompetent / unintelligent / competent / submissive
  # all sit close together).
  label_pos <- c(
    friendly      = 2L, unfriendly    = 2L, trust         = 2L,
    untrust       = 4L, intelligent   = 2L, dominant      = 4L,
    submissive    = 2L, competent     = 4L,
    incompetent   = 2L, unintelligent = 1L
  )[d$trait]
  graphics::text(d$icc_3_k, d$group_z, labels = d$trait,
                 pos = label_pos, cex = 0.75, col = "grey20",
                 offset = 0.5)

  # Annotation: r and 95% CI rounded to 2dp; p to 3dp, with the
  # "< 0.001" convention when the rounded value would be 0.000.
  p_str <- if (cor_stats[["p_value"]] < 0.001) "< 0.001" else
    sprintf("= %.3f", cor_stats[["p_value"]])
  ann <- sprintf(
    "Pearson r = %.2f, 95%% CI [%.2f, %.2f]\n(t(%d) = %.2f, p %s, n = 10 traits)",
    cor_stats[["estimate"]], cor_stats[["ci_lo"]], cor_stats[["ci_hi"]],
    as.integer(cor_stats[["df"]]),
    cor_stats[["t"]], p_str
  )
  graphics::legend("topleft", legend = ann, bty = "n",
                   text.col = "grey20", cex = 0.78,
                   inset = c(0.02, 0.01))
}

scatter_icc_vs_groupz(icc_vs_groupz$data,
                      icc_vs_groupz$pred,
                      icc_vs_groupz$cor,
                      file.path(fig_dir, "icc_vs_groupz.png"))

message("Done.")
message(sprintf(
  "ICC(3,k) vs group-mean z: r = %.3f, 95%% CI [%.3f, %.3f], p = %.2g",
  icc_vs_groupz$cor[["estimate"]],
  icc_vs_groupz$cor[["ci_lo"]],
  icc_vs_groupz$cor[["ci_hi"]],
  icc_vs_groupz$cor[["p_value"]]
))
