## Face-mask demonstration figures shown in vignette section 4.5.
##
## Inputs (under temp/base_faces/):
##   - FMNES.jpg                                 (RGB, 512x512)
##   - premade mask images/FMNES_mask_fullface_256.jpg (256x256)
##   - base_face_artficialperson.jpg             (grayscale, 512x512)
##
## Outputs (under vignettes/figures/face_masks/):
##   - fmnes_raw.png, fmnes_masked.png
##   - artificial_full.png, artificial_eyes.png,
##     artificial_nose.png, artificial_mouth.png,
##     artificial_upper_face.png, artificial_lower_face.png

suppressPackageStartupMessages({
  library(devtools)
})
devtools::load_all(".")

src_dir <- "temp/base_faces"
out_dir <- "vignettes/figures/face_masks"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

img_dims <- c(256L, 256L)

resize_nn <- function(mat, target = 256L) {
  if (identical(dim(mat), c(target, target))) return(mat)
  r <- nrow(mat); c <- ncol(mat)
  ri <- round(seq(1, r, length.out = target))
  ci <- round(seq(1, c, length.out = target))
  mat[ri, ci]
}

read_gray_resized <- function(path, target = 256L) {
  img <- jpeg::readJPEG(path)
  if (length(dim(img)) == 3L) {
    img <- 0.2126 * img[, , 1L] +
             0.7152 * img[, , 2L] +
             0.0722 * img[, , 3L]
  }
  img <- pmin(pmax(img, 0), 1)
  resize_nn(img, target)
}

# ---- read FMNES base + mask, resize to 256 ----------------------

message("Reading FMNES base + mask ...")
fmnes_raw  <- read_gray_resized(file.path(src_dir, "FMNES.jpg"), 256L)
fmnes_mask <- read_gray_resized(
  file.path(src_dir, "premade mask images",
            "FMNES_mask_fullface_256.jpg"),
  256L
)
fmnes_mask_logical <- fmnes_mask > 0.5

# Save raw FMNES as PNG so plot_ci_overlay can use it as a base.
fmnes_raw_path <- file.path(out_dir, "fmnes_raw_base.png")
png::writePNG(fmnes_raw, fmnes_raw_path)

# ---- FMNES showcase ---------------------------------------------

# Panel A: raw FMNES (no mask).
grDevices::png(file.path(out_dir, "fmnes_raw.png"),
               width = 600, height = 600, res = 100)
op <- graphics::par(mar = c(0.5, 0.5, 2, 0.5))
graphics::image(
  t(fmnes_raw)[, nrow(fmnes_raw):1],
  col      = grDevices::gray.colors(256, start = 0, end = 1),
  asp      = 1,
  axes     = FALSE,
  main     = "FMNES base face (raw)"
)
graphics::par(op)
grDevices::dev.off()

# Panel B: only the inside-mask region of the face is visible;
# outside is dimmed to mid-grey so the reader sees what the
# analysis "looks at".
fmnes_isolated <- fmnes_raw
fmnes_isolated[!fmnes_mask_logical] <- 0.85   # outside-mask -> light grey
grDevices::png(file.path(out_dir, "fmnes_masked.png"),
               width = 600, height = 600, res = 100)
op <- graphics::par(mar = c(0.5, 0.5, 2, 0.5))
graphics::image(
  t(fmnes_isolated)[, nrow(fmnes_isolated):1],
  col      = grDevices::gray.colors(256, start = 0, end = 1),
  asp      = 1,
  axes     = FALSE,
  main     = "Same face, restricted to a premade mask"
)
graphics::par(op)
grDevices::dev.off()

# ---- artificial-person base + 6 region masks --------------------

message("Reading artificial-person base + rendering 6 regions ...")
artif_raw <- read_gray_resized(
  file.path(src_dir, "base_face_artficialperson.jpg"),
  256L
)
artif_path <- file.path(out_dir, "artif_base.png")
png::writePNG(artif_raw, artif_path)

regions <- c("full", "eyes", "nose",
             "mouth", "upper_face", "lower_face")

for (region in regions) {
  m <- make_face_mask(img_dims, region = region)
  out <- file.path(out_dir,
                   sprintf("artificial_%s.png", region))
  grDevices::png(out, width = 600, height = 600, res = 100)
  plot_face_mask(
    m,
    img_dims   = img_dims,
    base_image = artif_path,
    main       = paste0("region = \"", region, "\"")
  )
  grDevices::dev.off()
}

message("Face-mask demo figures complete: ", out_dir)
