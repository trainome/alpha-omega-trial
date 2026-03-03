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
  "#000000"
)


alpha_theme <- function(yaxis = TRUE) {
  if (yaxis) {
    fun <- theme_classic() +
      theme(strip.background = element_blank())
  }
  if (!yaxis) {
    fun <- theme_classic() +
      theme(
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank()
      )
  }

  return(fun)
}
