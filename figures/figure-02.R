### Study design figure ###################################################

# Training data analysis

traindat <- AOData::training |>

  filter(!(participant %in% c(593, 613))) |>

  summarise(
    .by = c(participant, session),
    volume = sum(volume, na.rm = TRUE)
  ) |>
  inner_join(AOData::idkeys) |>

  filter(session != 21) |>
  mutate(baseline = if_else(session == 1, volume, NA)) |>
  mutate(.by = participant, baseline = sum(baseline, na.rm = TRUE)) |>
  mutate(rel_training = 100 * ((volume / baseline) - 1))


training_volume_fig <- traindat |>
  summarise(
    .by = c(session, treat),
    m = mean(rel_training),
    s = sd(rel_training)
  ) |>
  ggplot(aes(session, m)) +
  geom_smooth(color = colscale[1]) +
  alpha_theme() +
  labs(x = "Session", y = "Change in training volume<br>(% from session 1)") +
  scale_y_continuous(limits = c(-5, 50), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1, 20), expand = c(0, 0)) +
  theme(axis.title.y = element_markdown(angle = 0))


# Map study week to x position (week 0 at 0.14, week 20 at 0.9)
wk <- function(w) 0.14 + 0.76 * w / 20

# Left edge of the enrollment and allocation box (right edge at week 0)
x_enrol <- 0.02

# Top of the canvas (y units)
y_max <- 1.19

# Height of the horizontal dashed line connecting the measurement box to the
# measurement time-points
y_bend <- 0.86

# Study periods: enrollment and allocation, supplementation, lead-in and
# combined period. Each box ends at the week given in `end`.
periods <- tibble(
  xmin = c(x_enrol, wk(c(0, 7, 10))),
  xmax = wk(c(0, 7, 10, 20)),
  fill = c("white", "#E6E6E6", "#BAE3D8", "#4EBC9E"),
  label = c(
    "Enrollment and\nallocation",
    "Supplementation period,\nKrill oil vs. Placebo\n(7 weeks)",
    "Lead-in\n(3 weeks)",
    "Combined period,\nSupplementation (Krill oil vs. Placebo)\nand Resistance training (10 weeks)"
  ),
  end = c("Week 0", "Week 7", "Week 10", "Week 20")
)

# Measurement time-points (weeks)
meas_wk <- c(0, 7, 20)

# Measurement box, aligned with the period boxes
x_meas <- c(x_enrol, wk(20))

# Measurements grouped by collection method, placed in one row inside the
# measurement box (x: left edge of each column relative to the box, y: top of
# each group). The x positions give ~6 mm between columns at 7 pt text and
# fig_w = 180 mm.
measurements <- tribble(
  ~x,
  ~y,
  ~method,
  ~outcomes,
  0.024,
  1.14,
  "Blood samples (fasting)",
  "Lipids (TG, HDL, LDL),\nOmega-3 index, fatty acid profile\nAdiponectin",
  0.258,
  1.14,
  "OGTT (0-120 min)",
  "Glucose, C-peptide (0, 60 min)\nHOMA2-IR",
  0.483,
  1.14,
  "DXA",
  "Relative Fat mass",
  0.626,
  1.14,
  "Anthropometry and blood pressure",
  "Body mass, BMI\nWaist circumference\nSystolic and diastolic blood pressure",
  0.483,
  1.06,
  "Diet",
  "Food frequency\nquestionaire"
) |>
  mutate(x = x_meas[1] + x)

study_design <- data.frame(x = c(0, 1), y = c(0, 1)) |>

  ggplot(aes(x, y)) +

  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0)) +

  # No expansion, data coordinates match the cowplot canvas (0-1, y scaled
  # by y_max)
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +

  # Enrollment and allocation, study periods
  annotate(
    "rect",
    xmin = periods$xmin,
    xmax = periods$xmax,
    ymin = 0.53,
    ymax = 0.7,
    fill = periods$fill,
    alpha = 0.2,
    color = "black",
    linewidth = 0.3
  ) +
  annotate(
    "text",
    x = (periods$xmin + periods$xmax) / 2,
    y = 0.615,
    label = periods$label,
    size = 7 / .pt
  ) +

  # Week at the end of each period box
  annotate(
    "text",
    x = periods$xmax,
    y = 0.505,
    label = periods$end,
    size = 7 / .pt
  ) +

  # Measurements
  annotate(
    "rect",
    xmin = x_meas[1],
    xmax = x_meas[2],
    ymin = 0.95,
    ymax = 1.16,
    fill = "gray99",
    color = "black",
    linewidth = 0.3
  ) +
  annotate(
    "text",
    x = measurements$x,
    y = measurements$y,
    label = measurements$method,
    fontface = "bold",
    hjust = 0,
    vjust = 1,
    size = 7 / .pt
  ) +
  annotate(
    "text",
    x = measurements$x,
    y = measurements$y - 0.04,
    label = measurements$outcomes,
    hjust = 0,
    vjust = 1,
    lineheight = 1,
    size = 7 / .pt
  ) +

  # Measurement time-points, short names as in other figures. A dashed
  # vertical line from the centre of the measurement box down to y_bend, a
  # horizontal line spanning the time-points, then vertical arrows down to
  # just above the period boxes.
  annotate(
    "segment",
    x = mean(x_meas),
    xend = mean(x_meas),
    y = 0.95,
    yend = y_bend,
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  annotate(
    "segment",
    x = wk(min(meas_wk)),
    xend = wk(max(meas_wk)),
    y = y_bend,
    yend = y_bend,
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  annotate(
    "segment",
    x = wk(meas_wk),
    xend = wk(meas_wk),
    y = y_bend,
    yend = 0.73,
    linewidth = 0.3,
    linetype = "dashed",
    arrow = arrow(length = unit(1.5, "mm"), type = "closed")
  ) +
  annotate(
    "text",
    x = wk(meas_wk) + 0.01,
    y = (y_bend + 0.73) / 2,
    label = c("Baseline", "Supplementation", "Combined"),
    hjust = 0,
    size = 7 / .pt
  )

# Insert training volume figure under the combined period. The inset is
# widened by the space taken up left (y-axis) and right of its panel so that
# session 1 aligns with week 10 and session 20 with week 20. `fig_w` must match
# the width (mm) of the saved figure.
fig_w <- 180
train_grob <- ggplotGrob(training_volume_fig)
panel_col <- unique(train_grob$layout$l[train_grob$layout$name == "panel"])
panel_side_w <- function(cols) {
  sum(grid::convertWidth(train_grob$widths[cols], "mm", valueOnly = TRUE)) /
    fig_w
}
left_w <- panel_side_w(seq_len(panel_col - 1))
right_w <- panel_side_w((panel_col + 1):length(train_grob$widths))

study_design_fig <- cowplot::ggdraw(study_design) +
  cowplot::draw_plot(
    training_volume_fig,
    x = wk(10) - left_w,
    y = 0.04 / y_max,
    width = wk(20) - wk(10) + left_w + right_w,
    height = 0.42 / y_max
  )

# Missing data figure
# At each time-point a percentage of the complete data set is missing.
# This is a reason for using Bayesian imputation. We could add the percentage of
# observed data in the figure here...

# Combined figure: CONSORT (A) above study design (B). Panel heights (mm)
# keep 1.6 mm per canvas unit in the CONSORT and 80 mm per y unit in the
# study design.

design_h <- 95

figure_02 <- study_design_fig


ggsave(
  figure_02,
  filename = "figures/figure-2.pdf",
  width = fig_w,
  height = design_h,
  unit = "mm"
)

ggsave(
  figure_02,
  filename = "figures/figure-2.svg",
  width = fig_w,
  height = design_h,
  unit = "mm"
)


saveRDS(figure_02, "figures/figure-02.RDS")
