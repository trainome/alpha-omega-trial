# Enrollment figure and study design
#
# CONSORT flow diagram drawn with plain ggplot2 (geom_rect, geom_text,
# geom_segment). Boxes are defined in a data frame with centre coordinates
# on a 0-100 x 0-100 canvas; box heights are computed from the number of text
# lines. Fill in the counts in `n` (NA is printed as "XX") and adjust the
# box positions (x, y) if labels grow.
#
# Text size and line heights are tuned for the ggsave() size at the bottom
# (180 x 200 mm); if the output size changes, adjust `line_h` and `txt_size`.

# Source figure options, packages
source("./R/fig-opt-fun.R")
source("./R/packages.R")


# Checking number of participants in the data
# Number of randomized individuals. This collects all participants for
# which we have baseline data, these are regarded as "randomized"
part <- AOData::idkeys

## Descriptives data
baseline_dat <- AOData::bodycomp |>
  select(participant, time, lean.whole, BMC.whole, fat.whole) |>
  mutate(weight = lean.whole + BMC.whole + fat.whole) |>
  inner_join(part) |>
  filter(time == "T1") |>
  ## Calculate height based on weight (DXA) and BMI
  mutate(height = 100 * sqrt((weight / 1000) / bmi)) |>
  ## Calculate relative lean and fat
  mutate(
    lean.rel = 100 * (lean.whole / weight),
    fat.rel = 100 * (fat.whole / weight),
    # Weight in kg
    weight = weight / 1000
  ) |>
  select(participant, treat, completed, sex, age, height, weight, bmi) |>
  pivot_longer(cols = age:bmi)


## Included participants
totaln_rando <- baseline_dat |>
  distinct(participant) |>
  nrow()


# Allocated to placebo / N-3
allocation <- baseline_dat |>
  distinct(treat, participant) |>

  summarise(.by = treat, n = n()) |>
  with(setNames(n, treat))

# Lost to t2
t2_n <- AOData::idkeys |>
  filter(!is.na(round)) |>
  filter(completed != "t1") |>
  summarise(.by = treat, n = n()) |>
  with(setNames(n, treat))

# Lost to t4
t4_n <- AOData::idkeys |>
  filter(!is.na(round)) |>
  filter(completed == "t4") |>
  summarise(.by = treat, n = n()) |>
  with(setNames(n, treat))

# Counts ######################################################################

n <- list(
  assessed = 319,
  excluded = 172,
  excl_criteria = 40,
  excl_declined = 132,
  excl_other = 0,
  randomised = totaln_rando,

  # Per arm: allocation, follow-up periods and analysis.
  # lost_*: discontinued the intervention and lost to follow-up.
  # analysed: participants included in the (imputation) models, i.e. all
  # randomised unless some are excluded; observed: completed final follow-up.
  pla_alloc = allocation["ctrl"],
  pla_received = allocation["ctrl"],
  pla_not_received = 0,
  pla_lost_t2 = allocation["ctrl"] - t2_n["ctrl"],
  pla_lost_t4 = t2_n["ctrl"] - t4_n["ctrl"],
  pla_analysed = allocation["ctrl"],
  pla_observed = t4_n["ctrl"],
  pla_excl_analysis = 0,

  n3_alloc = allocation["omega"],
  n3_received = allocation["omega"],
  n3_not_received = 0,
  n3_lost_t2 = allocation["omega"] - t2_n["omega"],
  n3_lost_t4 = t2_n["omega"] - t4_n["omega"],
  n3_analysed = allocation["omega"],
  n3_observed = t4_n["omega"],
  n3_excl_analysis = 0
)

# Participants with imputed final follow-up data
n$pla_imputed <- n$pla_analysed - n$pla_observed
n$n3_imputed <- n$n3_analysed - n$n3_observed

# Reasons ######################################################################

# Reasons listed under each arm box, one string per line. Use character(0)
# to omit the reasons lines; add "\n" inside a string to wrap long reasons.
reasons <- list(
  pla_alloc = c(character(0)),
  pla_t2 = c(
    "Reasons not related to the study (n=3)",
    "Fear of COVID-19 (n=6)"
  ),
  pla_t4 = c(
    "Lockdown (n=9)",
    "Fear of COVID-19 (n=11)",
    "Low training adherence (n=2)"
  ),
  pla_analysis = c(character(0)),

  n3_alloc = c(character(0)),
  n3_t2 = c(
    "Reasons not related to the study (n=3)",
    "Fear of COVID-19 (n=7)",
    "Not enough time (n=2)"
  ),
  n3_t4 = c(
    "Lockdown (n=7)",
    "Fear of COVID-19 (n=5)",
    "Low training adherence (n=1)"
  ),
  n3_analysis = c(character(0))
)

# Format counts, NA -> placeholder
f <- function(x) if (is.na(x)) "XX" else format(x, big.mark = ",")


# Box labels ##################################################################

arm_boxes <- function(arm, prefix) {
  g <- function(x) f(n[[paste0(prefix, "_", x)]])

  # Indented reasons lines (empty string if no reasons are given)
  r <- function(x) {
    rs <- reasons[[paste0(prefix, "_", x)]]
    if (length(rs) == 0) "" else paste0("\n    ", rs, collapse = "")
  }

  tibble(
    id = paste0(prefix, c("_alloc", "_t2", "_t4", "_analysis")),
    label = c(
      paste0(
        "Allocated to ",
        arm,
        " (n = ",
        g("alloc"),
        ")\n",
        "  Received allocated intervention (n = ",
        g("received"),
        ")\n",
        "  Did not receive allocated intervention (n = ",
        g("not_received"),
        ")",
        r("alloc")
      ),
      paste0(
        "Supplementation period (Week 1 to week 7)\n",
        "  Discontinued and lost to follow-up (n = ",
        g("lost_t2"),
        ")",
        r("t2")
      ),
      paste0(
        "Supplementation + training period (Week 8 to week 23)\n",
        "  Discontinued and lost to follow-up (n = ",
        g("lost_t4"),
        ")",
        r("t4")
      ),
      paste0(
        "Analysed (n = ",
        g("analysed"),
        ")\n",
        "  Observed at final follow-up (n = ",
        g("observed"),
        ")\n",
        "  Missing follow-up data imputed (n = ",
        g("imputed"),
        ")\n",
        "  Excluded from analysis (n = ",
        g("excl_analysis"),
        ")",
        r("analysis")
      )
    )
  )
}


# Layout ######################################################################

line_h <- 1.6 # height of one text line (canvas units)
pad <- 1.5 # vertical padding inside boxes
txt_size <- 2.4 # text size (mm)
y_lim <- c(4, 89.5) # vertical extent of the canvas (crops empty space)

x_main <- 54.5 # centre of the main (enrollment) column
x_arm <- c(pla = 34, n3 = 75) # centres of the two arm columns
w_arm <- 38 # width of arm boxes

boxes <- bind_rows(
  tibble(
    id = c("assessed", "excluded", "randomised"),
    label = c(
      paste0("Assessed for eligibility (n = ", f(n$assessed), ")"),
      paste0(
        "Excluded (n = ",
        f(n$excluded),
        ")\n",
        "  Not meeting inclusion criteria (n = ",
        f(n$excl_criteria),
        ")\n",
        "  Declined to participate (n = ",
        f(n$excl_declined),
        ")\n",
        "  Other reasons (n = ",
        f(n$excl_other),
        ")"
      ),
      paste0("Randomised (n = ", f(n$randomised), ")")
    ),
    x = c(x_main, 81, x_main),
    y = c(87, 79, 71),
    w = c(36, 32, 30)
  ),
  arm_boxes("Placebo", "pla") |>
    mutate(x = x_arm[["pla"]], y = c(57, 41, 25, 9), w = w_arm),
  arm_boxes("N-3", "n3") |>
    mutate(x = x_arm[["n3"]], y = c(57, 41, 25, 9), w = w_arm)
) |>
  mutate(
    n_lines = str_count(label, "\n") + 1,
    h = n_lines * line_h + pad,
    xmin = x - w / 2,
    xmax = x + w / 2,
    ymin = y - h / 2,
    ymax = y + h / 2,
    # Single-line boxes are centred, multi-line boxes left-aligned
    hjust = if_else(n_lines == 1, 0.5, 0),
    x_txt = if_else(n_lines == 1, x, xmin + 1.5)
  )

# Box lookup helper
b <- function(id, col) boxes[[col]][boxes$id == id]


# Arrows ######################################################################

# Split point between randomisation and allocation
y_split <- (b("randomised", "ymin") + b("pla_alloc", "ymax")) / 2

# Connecting lines without arrow heads
lines <- tribble(
  ~x,
  ~y,
  ~xend,
  ~yend,
  # Randomised -> split
  x_main,
  b("randomised", "ymin"),
  x_main,
  y_split,
  # Horizontal split line
  x_arm[["pla"]],
  y_split,
  x_arm[["n3"]],
  y_split
)

# Vertical arrows within each arm
arm_chain <- function(prefix) {
  ids <- paste0(prefix, c("_alloc", "_t2", "_t4", "_analysis"))
  tibble(
    x = b(ids[1], "x"),
    y = sapply(ids[-4], b, col = "ymin"),
    xend = b(ids[1], "x"),
    yend = sapply(ids[-1], b, col = "ymax")
  )
}

arrows <- bind_rows(
  # Assessed -> randomised
  tibble(
    x = x_main,
    y = b("assessed", "ymin"),
    xend = x_main,
    yend = b("randomised", "ymax")
  ),
  # Main column -> excluded
  tibble(
    x = x_main,
    y = b("excluded", "y"),
    xend = b("excluded", "xmin"),
    yend = b("excluded", "y")
  ),
  # Split -> allocation boxes
  tibble(
    x = unname(x_arm),
    y = y_split,
    xend = unname(x_arm),
    yend = b("pla_alloc", "ymax")
  ),
  arm_chain("pla"),
  arm_chain("n3")
)


# Stage labels (left margin) ##################################################

stages <- tribble(
  ~stage,
  ~ymin,
  ~ymax,
  "Enrollment",
  b("randomised", "ymin"),
  b("assessed", "ymax"),
  "Allocation",
  b("pla_alloc", "ymin"),
  b("pla_alloc", "ymax"),
  "Follow-up",
  b("pla_t4", "ymin"),
  b("pla_t2", "ymax"),
  "Analysis",
  b("pla_analysis", "ymin"),
  b("pla_analysis", "ymax")
) |>
  mutate(xmin = 1, xmax = 12)


# Plot ########################################################################

consort <- ggplot() +
  geom_rect(
    data = stages,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "gray95",
    color = "gray80",
    linewidth = 0.3
  ) +
  geom_text(
    data = stages,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = stage),
    size = txt_size,
    fontface = "bold"
  ) +
  geom_segment(
    data = lines,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.3
  ) +
  geom_segment(
    data = arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.3,
    arrow = arrow(length = unit(1.5, "mm"), type = "closed")
  ) +
  geom_rect(
    data = boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    data = boxes,
    aes(x = x_txt, y = y, label = label, hjust = hjust),
    size = txt_size,
    lineheight = 1
  ) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  theme_void()


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
  labs(x = "Session", y = "Training volume<br>(% of session 1)") +
  scale_y_continuous(limits = c(-5, 50), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1, 20), expand = c(0, 0)) +
  theme(axis.title.y = element_markdown(angle = 0))


# Map study week to x position (week 0 at 0.14, week 20 at 0.9)
wk <- function(w) 0.14 + 0.76 * w / 20

# Left edge of the enrollment and allocation box (right edge at week 0)
x_enrol <- 0.02

# Top of the canvas (y units; the time line is at 0.1)
y_max <- 1.19

# Height where the measurement arrows bend from angled to vertical
y_bend <- 0.86

# Measurements grouped by collection method, placed in one row inside the
# measurement box (x: left edge of each column, y: top of each group). The
# x positions give ~6 mm between columns at 7 pt text and fig_w = 180 mm.
measurements <- tribble(
  ~x,
  ~y,
  ~method,
  ~outcomes,
  0.084,
  1.12,
  "Blood samples (fasting)",
  "Lipids (TG, HDL, LDL),\nOmega-3 index, fatty acid profile\nAdiponectin",
  0.318,
  1.12,
  "OGTT (0-120 min)",
  "Glucose, C-peptide (0, 60 min)\nHOMA2-IR",
  0.543,
  1.12,
  "DXA",
  "Relative Fat mass",
  0.686,
  1.12,
  "Anthropometry and blood pressure",
  "Body mass, BMI\nWaist circumference\nSystolic and diastolic blood pressure"
)

study_design <- data.frame(x = c(0, 1), y = c(0, 1)) |>

  ggplot(aes(x, y)) +

  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0)) +

  # No expansion, data coordinates match the cowplot canvas (0-1, y scaled
  # by y_max)
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +

  # Time scale
  annotate(
    "segment",
    x = wk(0),
    xend = wk(20),
    y = 0.1,
    yend = 0.1,
    linewidth = 0.3
  ) +
  annotate(
    "segment",
    x = wk(c(0, 8, 20)),
    xend = wk(c(0, 8, 20)),
    y = 0.08,
    linewidth = 0.3,
    yend = 0.1
  ) +

  annotate(
    "text",
    x = wk(c(0, 8, 20)),
    y = 0.05,
    label = c("Week 0", "Week 8", "Week 20"),
    size = 7 / .pt
  ) +

  # Enrollment and allocation, study periods
  annotate(
    "rect",
    xmin = c(x_enrol, wk(c(0, 8))),
    xmax = wk(c(0, 8, 20)),
    ymin = 0.53,
    ymax = 0.7,
    fill = c("white", colscale[c(2, 4)]),
    alpha = 0.2,
    color = "black",
    linewidth = 0.3
  ) +
  annotate(
    "text",
    x = c(x_enrol + 0.014, wk(c(2, 10))),
    y = 0.615,
    label = c(
      "Enrollment and\nallocation",
      "Supplementation period,\nN-3 vs. Placebo",
      "Combined period,\nSupplementation (N-3 vs. Placebo)\nand Resistance training"
    ),
    hjust = 0,
    size = 7 / .pt
  ) +

  # Measurements
  annotate(
    "rect",
    xmin = 0.05,
    xmax = 0.95,
    ymin = 0.95,
    ymax = 1.16,
    fill = "gray80",
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

  # Measurement time-points, short names as in other figures. Angled
  # segments from one point at the bottom of the measurement box to y_bend,
  # then vertical arrows down to just above the period boxes.
  annotate(
    "segment",
    x = wk(10),
    xend = wk(c(0, 8, 20)),
    y = 0.95,
    yend = y_bend,
    linewidth = 0.3
  ) +
  annotate(
    "segment",
    x = wk(c(0, 8, 20)),
    xend = wk(c(0, 8, 20)),
    y = y_bend,
    yend = 0.73,
    linewidth = 0.3,
    arrow = arrow(length = unit(1.5, "mm"), type = "closed")
  ) +
  annotate(
    "text",
    x = wk(c(0, 8, 20)) + 0.01,
    y = (y_bend + 0.73) / 2,
    label = c("Baseline", "Supplementation", "Combined"),
    hjust = 0,
    size = 7 / .pt
  )

# Insert training volume figure under the combined period. The inset is
# widened by the space taken up left (y-axis) and right of its panel so that
# session 1 aligns with week 8 and session 20 with week 20. `fig_w` must match
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
    x = wk(8) - left_w,
    y = 0.12 / y_max,
    width = wk(20) - wk(8) + left_w + right_w,
    height = 0.34 / y_max
  )

# Missing data figure
# At each time-point a percentage of the complete data set is missing.
# This is a reason for using Bayesian imputation. We could add the percentage of
# observed data in the figure here...

# Combined figure: CONSORT (A) above study design (B). Panel heights (mm)
# keep 1.6 mm per canvas unit in the CONSORT and 80 mm per y unit in the
# study design.
consort_h <- 1.6 * diff(y_lim)
design_h <- 80 * y_max

figure_01 <- cowplot::plot_grid(
  cowplot::plot_grid(NULL, consort, rel_widths = c(0.05, 1)),
  study_design_fig,
  ncol = 1,
  rel_heights = c(consort_h, design_h),
  labels = c("A", "B"),
  label_size = 10
)


ggsave(
  figure_01,
  filename = "figures/figure-1.pdf",
  width = fig_w,
  height = consort_h + design_h,
  unit = "mm"
)

ggsave(
  figure_01,
  filename = "figures/figure-1.svg",
  width = fig_w,
  height = consort_h + design_h,
  unit = "mm"
)


saveRDS(figure_01, "figures/figure-01.RDS")
