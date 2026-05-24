## Face-mask demonstration figures shown in vignette section 4.5.
##
## Inputs (under temp/base_faces/):
##   - FMNES.jpg                                 (RGB, 512x512)
##   - premade mask images/FMNES_mask_fullface_256.jpg (256x256)
##   - base_face_artficialperson.jpg             (grayscale, 512x512)
##
## Outputs (under vignettes/figures/face_masks/):
##   - fmnes_raw.png, fmnes_masked.png
##   - artificial_full.png, artificial_upper_face.png,
##     artificial_lower_face.png, artificial_nose.png,
##     artificial_mouth.png, artificial_eyes.png,
##     artificial_left_eye.png, artificial_right_eye.png
##   - artificial_mouth_demo_default.png,
##     artificial_mouth_demo_shift_v.png,
##     artificial_mouth_demo_shift_vh.png  (shift_mask demo)
##   - artificial_left_eye_demo_default.png,
##     artificial_left_eye_demo_tuned.png  (region_bounds demo)

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

regions <- c("full", "upper_face", "lower_face", "nose",
             "mouth", "eyes", "left_eye", "right_eye")

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

# ---- shift_mask tuning demo on the ELLIPTICAL mouth region ----

message("Rendering shift_mask demo on the mouth region ...")
# Use the exported rcisignal::shift_mask (renamed args in v0.2.0).

mouth_default <- matrix(
  make_face_mask(img_dims, region = "mouth"),
  nrow = img_dims[1L], ncol = img_dims[2L]
)
# On this base face the actual mouth sits above the default mask
# position, so a real tune moves the mask up (positive vertical;
# see `?shift_mask` for the math / y-axis-up sign convention).
# Combine with a small rightward shift to demonstrate the two-axis
# form.
mouth_v  <- shift_mask(mouth_default, vertical = 20)
mouth_vh <- shift_mask(mouth_default, vertical = 20, horizontal = 8)

mouth_panels <- list(
  artificial_mouth_demo_default  = list(mask = mouth_default,
                                        main = "mouth: default"),
  artificial_mouth_demo_shift_v  = list(mask = mouth_v,
                                        main = "mouth: up 20 px"),
  artificial_mouth_demo_shift_vh = list(mask = mouth_vh,
                                        main = "mouth: up 20 + right 8 px")
)
for (nm in names(mouth_panels)) {
  out <- file.path(out_dir, sprintf("%s.png", nm))
  grDevices::png(out, width = 600, height = 600, res = 100)
  plot_face_mask(
    mouth_panels[[nm]]$mask,
    img_dims   = img_dims,
    base_image = artif_path,
    main       = mouth_panels[[nm]]$main
  )
  grDevices::dev.off()
}

# ---- region_bounds tuning demo on the left_eye rectangle ------

message("Rendering region_bounds demo on the left_eye region ...")
le_default <- make_face_mask(img_dims, region = "left_eye")
le_tuned   <- make_face_mask(
  img_dims, region = "left_eye",
  region_bounds = c(0.40, 0.50, 0.24, 0.44)
)

le_panels <- list(
  artificial_left_eye_demo_default = list(mask = le_default,
                                          main = "left_eye: default"),
  artificial_left_eye_demo_tuned   = list(mask = le_tuned,
                                          main = "left_eye: tuned bounds")
)
for (nm in names(le_panels)) {
  out <- file.path(out_dir, sprintf("%s.png", nm))
  grDevices::png(out, width = 600, height = 600, res = 100)
  plot_face_mask(
    le_panels[[nm]]$mask,
    img_dims   = img_dims,
    base_image = artif_path,
    main       = le_panels[[nm]]$main
  )
  grDevices::dev.off()
}

message("Face-mask demo figures complete: ", out_dir)
