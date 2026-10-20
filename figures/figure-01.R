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

x_main <- 50 # centre of the main (enrollment) column
x_arm <- c(pla = 24, n3 = 76) # centres of the two arm columns
w_arm <- 46 # width of arm boxes

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
    x = c(x_main, 82, x_main),
    y = c(87, 79, 71),
    w = c(36, 32, 30)
  ),
  arm_boxes("Placebo", "pla") |>
    mutate(x = x_arm[["pla"]], y = c(57, 41, 25, 9), w = w_arm),
  arm_boxes("Krill-oil", "n3") |>
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


# Stage labels ################################################################

# Enrollment sits level with "Assessed for eligibility", its left edge aligned
# with the placebo allocation box. The other stages are centred between the
# arms, midway in the gap above the first box of each step (using the taller
# box of the two arms).
stage_h <- line_h + pad # same height as a single-line box
w_stage <- 14

# Top of the highest box in a step across both arms
step_top <- function(step) max(sapply(paste0(c("pla", "n3"), step), b, col = "ymax"))
step_bot <- function(step) min(sapply(paste0(c("pla", "n3"), step), b, col = "ymin"))

stages <- tibble(
  stage = c("Enrollment", "Allocation", "Follow-up", "Analysis"),
  x = c(b("pla_alloc", "xmin") + w_stage / 2, x_main, x_main, x_main),
  y = c(
    b("assessed", "y"),
    (y_split + step_top("_alloc")) / 2,
    (step_bot("_alloc") + step_top("_t2")) / 2,
    (step_bot("_t4") + step_top("_analysis")) / 2
  ),
  w = w_stage
) |>
  mutate(
    xmin = x - w / 2,
    xmax = x + w / 2,
    ymin = y - stage_h / 2,
    ymax = y + stage_h / 2
  )

# Enrollment label -> assessed box
lines <- bind_rows(
  lines,
  tibble(
    x = stages$xmax[stages$stage == "Enrollment"],
    y = b("assessed", "y"),
    xend = b("assessed", "xmin"),
    yend = b("assessed", "y")
  )
)


# Plot ########################################################################

consort <- ggplot() +
  geom_rect(
    data = stages,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#d6e6f5",
    color = NA
  ) +
  geom_text(
    data = stages,
    aes(x = x, y = y, label = stage),
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


figure_01 <- consort

# Settings #################################################

fig_w <- 180 # Figure width
consort_h <- 140 # figure height
# Compile figure ###########################################

ggsave(
  figure_01,
  filename = "figures/figure-1.pdf",
  width = fig_w,
  height = consort_h,
  unit = "mm"
)

ggsave(
  figure_01,
  filename = "figures/figure-1.svg",
  width = fig_w,
  height = consort_h,
  unit = "mm"
)


saveRDS(figure_01, "figures/figure-01.RDS")
