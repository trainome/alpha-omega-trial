## Plotting functions and settings #######################

# cowplot measures text on a separate null device. With the default device,
# ggtext labels lose spaces when figures are drawn with ragg (as in the book);
# measure with ragg as well.
cowplot::set_null_device("agg")

# Color scale
colscale <- c(
  "#F0E442",
  "#D55E00", # Krill-oil
  "#0072B2", # Supplementation only
  "#999999", # Placebo
  "#0072B2", #
  "#CC79A7", # Combined
  "#D55E00",

  "darkmagenta",
  "coral2",
  "forestgreen",
  "deeppink2",
  "#000000"
)


alpha_theme <- function(yaxis = TRUE, text_base_size = 8) {
  if (yaxis) {
    fun <- theme_classic(base_size = text_base_size) +
      theme(
        strip.background = element_rect(fill = "gray95", color = "white"),
        strip.text = element_text(hjust = 0),
        axis.title.x = ggtext::element_markdown(size = text_base_size),
        axis.title.y = ggtext::element_markdown(size = text_base_size)
      )
  }
  if (!yaxis) {
    fun <- theme_classic(base_size = text_base_size) +
      theme(
        strip.background = element_rect(fill = "gray95", color = "white"),
        strip.text = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = ggtext::element_markdown(size = text_base_size)
      )
  }

  return(fun)
}
