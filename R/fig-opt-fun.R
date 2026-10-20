## Plotting functions and settings #######################

# Color scale
colscale <- c(
  "#009E73",
  "#D55E00",
  "#E69F00",
  "#56B4E9",
  "#F0E442",
  "#0072B2",
  "#CC79A7",
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
