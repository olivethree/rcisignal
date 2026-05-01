# Build the rcisignal hex logo.
#
# Run from the package root:
#   source("data-raw/make_logo.R")
# Output: man/figures/logo.png (and logo.svg) plus two alternate palettes.

pkgs <- c("hexSticker", "showtext", "ggplot2", "sysfonts")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install first: install.packages(c(",
       paste0("\"", missing, "\"", collapse = ", "), "))")
}

library(hexSticker)
library(ggplot2)
library(showtext)

sysfonts::font_add_google("JetBrains Mono", "jbm")
showtext_auto()

palettes <- list(
  cyan = list(
    bg     = "#0B1220",
    border = "#38BDF8",
    wave   = "#7DD3FC",
    ghost  = "#1E293B",
    text   = "#E0F2FE"
  ),
  ember = list(
    bg     = "#1A0A06",
    border = "#B22222",
    wave   = "#EEB422",
    ghost  = "#3B1A0E",
    text   = "#FFE9B0"
  ),
  cyberpunk = list(
    bg     = "#0A0014",
    border = "#FF006E",
    wave   = "#00F5FF",
    ghost  = "#3A0A4D",
    text   = "#F5D5FF"
  ),
  thermal = list(
    bg     = "#1B0C3F",
    border = "#F2580C",
    wave   = "#FCFFA4",
    ghost  = "#7C1D6F",
    text   = "#FED976"
  ),
  aurora = list(
    bg     = "#0A1E3F",
    border = "#14B8A6",
    wave   = "#C084FC",
    ghost  = "#1E3A8A",
    text   = "#E0E7FF"
  ),
  paper = list(
    bg     = "#F5F1E8",
    border = "#1F2937",
    wave   = "#DC2626",
    ghost  = "#D1D5DB",
    text   = "#1F2937"
  ),
  matrix = list(
    bg     = "#000000",
    border = "#00FF41",
    wave   = "#39FF14",
    ghost  = "#003B00",
    text   = "#39FF14"
  ),
  bubblegum = list(
    bg     = "#3D0066",
    border = "#FF4FB7",
    wave   = "#7FFFD4",
    ghost  = "#6B2B8C",
    text   = "#FFE4F1"
  ),
  galaxy = list(
    bg     = "#0B0033",
    border = "#FFD700",
    wave   = "#F8F8FF",
    ghost  = "#5B2A86",
    text   = "#FFD700"
  )
)

x <- seq(-pi, pi, length.out = 600)
envelope <- exp(-0.18 * x^2)
df_main  <- data.frame(x = x, y = sin(3 * x) * envelope)
df_ghost <- data.frame(x = x, y = sin(7 * x) * envelope * 0.35)

build_subplot <- function(p) {
  ggplot() +
    geom_hline(yintercept = 0, colour = p$ghost, linewidth = 0.4) +
    geom_line(data = df_ghost, aes(x, y),
              colour = p$ghost, linewidth = 0.8) +
    geom_line(data = df_main, aes(x, y),
              colour = p$wave, linewidth = 1.6, lineend = "round") +
    coord_cartesian(ylim = c(-1.2, 1.2)) +
    theme_void() +
    theme_transparent()
}

# showtext + the default SVG device disagree on text metrics, so the same
# `p_size` value renders differently in PNG vs SVG. Tune SVG sizes separately.
is_svg <- function(f) grepl("\\.svg$", f)

write_sticker <- function(p, filename) {
  args <- list(
    subplot   = build_subplot(p),
    package   = "rcisignal",
    p_size    = if (is_svg(filename)) 6   else 30,
    p_y       = 1.45,
    p_color   = p$text,
    p_family  = "jbm",
    s_x       = 1, s_y = 0.85,
    s_width   = 1.5, s_height = 0.9,
    h_fill    = p$bg,
    h_color   = p$border,
    h_size    = 1.4,
    url       = "",
    u_size    = if (is_svg(filename)) 1.5 else 4.2,
    u_color   = p$text,
    u_family  = "jbm",
    white_around_sticker = TRUE,
    filename  = filename
  )
  if (!is_svg(filename)) args$dpi <- 600
  do.call(sticker, args)
}

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

# `aurora` is the canonical logo; the others land alongside as alternates.
write_sticker(palettes$aurora,    "man/figures/logo.png")
write_sticker(palettes$aurora,    "man/figures/logo.svg")
write_sticker(palettes$cyan,      "man/figures/logo-cyan.png")
write_sticker(palettes$cyan,      "man/figures/logo-cyan.svg")
write_sticker(palettes$ember,     "man/figures/logo-ember.png")
write_sticker(palettes$ember,     "man/figures/logo-ember.svg")
write_sticker(palettes$cyberpunk, "man/figures/logo-cyberpunk.png")
write_sticker(palettes$cyberpunk, "man/figures/logo-cyberpunk.svg")
write_sticker(palettes$thermal,   "man/figures/logo-thermal.png")
write_sticker(palettes$thermal,   "man/figures/logo-thermal.svg")
write_sticker(palettes$aurora,    "man/figures/logo-aurora.png")
write_sticker(palettes$aurora,    "man/figures/logo-aurora.svg")
write_sticker(palettes$paper,     "man/figures/logo-paper.png")
write_sticker(palettes$paper,     "man/figures/logo-paper.svg")
write_sticker(palettes$matrix,    "man/figures/logo-matrix.png")
write_sticker(palettes$matrix,    "man/figures/logo-matrix.svg")
write_sticker(palettes$bubblegum, "man/figures/logo-bubblegum.png")
write_sticker(palettes$bubblegum, "man/figures/logo-bubblegum.svg")
write_sticker(palettes$galaxy,    "man/figures/logo-galaxy.png")
write_sticker(palettes$galaxy,    "man/figures/logo-galaxy.svg")

message("Wrote 18 files under man/figures/ (9 palettes x png+svg)")
