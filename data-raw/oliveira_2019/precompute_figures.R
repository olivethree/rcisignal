## Resume-only script: produce the figure outputs the main
## precompute.R does not generate.  Reads cached
## signal_matrices.rds + dissim_full.rds.

suppressPackageStartupMessages({
  library(devtools)
})
devtools::load_all(".")

data_dir <- Sys.getenv(
  "rcisignal_OLV_DATA",
  unset = "temp/rcicrdiagnostics/temp/oliveira_study"
)

fig_dir   <- "vignettes/figures/oliveira_2019"
cache_dir <- "vignettes/precomputed/oliveira_2019"

img_dims <- c(256L, 256L)
base_face_path <- file.path(fig_dir, "base_face.png")

sm          <- readRDS(file.path(cache_dir, "signal_matrices.rds"))
dissim_full <- readRDS(file.path(cache_dir, "dissim_full.rds"))

regions <- c("full", "eyes", "mouth", "upper_face")

# ---- 8. pairwise cluster-agreement maps -----------------------------

message("Generating pairwise cluster-agreement maps ...")

build_sig_mask <- function(ct) {
  sig_pos <- ct$clusters$cluster_id[
    ct$clusters$direction == "pos" & ct$clusters$significant
  ]
  sig_neg <- ct$clusters$cluster_id[
    ct$clusters$direction == "neg" & ct$clusters$significant
  ]
  as.vector((ct$pos_labels %in% sig_pos) |
              (ct$neg_labels %in% sig_neg))
}

pairwise <- list(
  list(name = "trust_vs_friendly",
       a = "trust", b = "friendly",
       title = "Trust minus Friendly (FWER-controlled clusters)"),
  list(name = "dominant_vs_competent",
       a = "dominant", b = "competent",
       title = "Dominant minus Competent (FWER-controlled clusters)")
)

oval_mask <- make_face_mask(img_dims, region = "full")

cluster_pairwise <- list()
for (p in pairwise) {
  message("  cluster test: ", p$name)
  ct <- rel_cluster_test(
    sm[[p$a]], sm[[p$b]],
    img_dims          = img_dims,
    cluster_threshold = 2.0,
    n_permutations    = 2000L,
    seed              = 1L,
    progress          = FALSE
  )
  cluster_pairwise[[p$name]] <- ct

  diff_signal <- rowMeans(sm[[p$a]]) - rowMeans(sm[[p$b]])

  # Descriptive map (no significance filter, restricted to oval).
  desc_path <- file.path(fig_dir,
                         sprintf("pairwise_descriptive_%s.png", p$name))
  desc_title <- sub(" \\(FWER-controlled clusters\\)$",
                    " (descriptive)", p$title)
  grDevices::png(desc_path, width = 600, height = 600, res = 100)
  plot_ci_overlay(
    diff_signal,
    base_image = base_face_path,
    img_dims   = img_dims,
    mask       = oval_mask,
    main       = desc_title
  )
  grDevices::dev.off()
  message("  wrote ", desc_path)

  # FWER-controlled map (significant-cluster mask).
  sig_mask <- build_sig_mask(ct)
  out_path <- file.path(fig_dir, sprintf("pairwise_%s.png", p$name))
  grDevices::png(out_path, width = 600, height = 600, res = 100)
  plot_ci_overlay(
    diff_signal,
    base_image = base_face_path,
    img_dims   = img_dims,
    mask       = sig_mask,
    main       = p$title
  )
  grDevices::dev.off()
  message("  wrote ", out_path)
}

saveRDS(cluster_pairwise,
        file.path(cache_dir, "cluster_pairwise.rds"))

# ---- 9. dissimilarity grid -----------------------------------------

message("Generating plot_dissimilarity_grid figure ...")
grDevices::png(file.path(fig_dir, "dissim_grid.png"),
               width = 900, height = 360, res = 100)
plot_dissimilarity_grid(
  "Trust vs Friendly"     = dissim_full[["Trust vs Friendly"]],
  "Competent vs Dominant" = dissim_full[["Competent vs Dominant"]],
  "Trust vs Dominant"     = dissim_full[["Trust vs Dominant"]]
)
grDevices::dev.off()

# ---- 10. face-mask reference figures -------------------------------

message("Generating face-mask reference figures ...")
for (region in regions) {
  m <- make_face_mask(img_dims, region = region)
  grDevices::png(
    file.path(fig_dir, sprintf("mask_%s.png", region)),
    width = 360, height = 360, res = 100
  )
  plot_face_mask(m, img_dims = img_dims,
                 base_image = base_face_path,
                 main = paste("region =", region))
  grDevices::dev.off()
}

message("Figure precompute complete.")
