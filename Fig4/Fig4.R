library(tidyverse)
library(ggtext)
library(ggbeeswarm)
library(cowplot)
library(here)

source(here("utils", "helper_functions.R"))

# import ----------------------------------------------------------------
# start at line 25

# long_df <- feathers |>
#   filter(date %in% c('20250626', '20250703', '20250717'))

# # keep any track > 10 seconds (15 FPS)
# long_frame_filter <- long_df |>
#   group_by(video, particle) |>
#   tally() |>
#   filter(n > 150)

# long_filtered <- long_frame_filter |>
#   select(-n) |>
#   left_join(long_df)

long_filtered <- read_rds(here("Fig4", "long_filtered.rds"))

# nest_cols <- c("date", "video", "particle")
# long_nested <- quick_nest(long_filtered, nest_cols)

# long_track_summary <- calculate_track_features_parallel(
#   long_nested,
#   fps = 15,
#   chunk_size = 500
# )

long_track_summary <- read_rds(here("Fig4", "long_track_summary.rds"))

long_tracks <- long_track_summary |>
  unnest(c(data)) |>
  mutate(
    camera = str_extract(video, "[0-9]{8}$"),
    tissue = case_when(
      camera %in% c("24568709", "24568744") ~ "intestine",
      camera %in% c("25128038", "25112214") ~ "liver"
    ),
    .before = particle
  ) |>
  # a bit of filtration
  filter(
    speed_mean > 5,
    net_displacement > 100,
    sd_x > 10 & sd_y > 10
  )

# binning

max_frames <- 4 * 60 * 60 * 15
interval <- 5 * 60 * 15 # 5 minutes in frames

(attrition_plot <- long_tracks |>
  ungroup() |>
  mutate(
    chunk = cut(
      frame,
      breaks = seq(0, max_frames, by = interval),
      include.lowest = TRUE,
      labels = FALSE
    ) *
      5
  ) |>
  group_by(video, tissue, chunk) |>
  tally() |>
  ggplot() +
  geom_smooth(
    aes(x = chunk, y = n, group = video, color = tissue),
    se = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(0, max_frames / 15 / 60, by = 30),
    labels = seq(0, max_frames / 15 / 60, by = 30)
  ) +
  scale_color_manual(
    values = c("intestine" = "indianred", "liver" = "steelblue"),
    labels = c("intestine" = "Intestine", "liver" = "Liver")
  ) +
  labs(
    x = "Time (minutes)",
    y = "Number of observations in 5 minute bin",
    color = "Tissue source"
  ) +
  theme_half_open() +
  theme(
    panel.grid = element_line(colour = "grey92", size = rel(0.5)),
    panel.grid.minor = element_blank(),
    axis.title = element_markdown(size = 9),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.title = element_markdown(size = 9),
    legend.text = element_markdown(size = 8),
  ) +
  NULL)

save_plot(
  here("Fig4", "Fig4.pdf"),
  attrition_plot,
  base_height = 4,
  base_width = 6
)

save_plot(
  here("Fig4", "Fig4.png"),
  attrition_plot,
  base_height = 4,
  base_width = 6,
  bg = 'white'
)

# long_chunked <- long_filtered |>
#   group_by(video, particle) |>
#   group_split() |>
#   map_dfr(~ split_trajectory(.x, frame_rate = 15, chunk_duration_sec = 5))

# lnest_cols <- c("date", "video", "particle", "subparticle")
# long_nested <- quick_nest(long_chunked, nest_cols)

# long_subtrack_summary <- calculate_track_features_parallel(
#   long_nested,
#   fps = 15,
#   chunk_size = 1000
# )

long_subtrack_summary <- read_rds(here('Fig4', 'long_subtrack_summary.rds')) |>
  mutate(
    date = str_extract(video, "2025[0-9]{4}"),
    camera = str_extract(video, "[0-9]{8}$"),
    # response (generic variable so we can use fit_model across experiments) == region of interest
    tissue = case_when(
      camera %in% c("24568709", "24568744") ~ "intestine",
      camera %in% c("25128038", "25112214") ~ "liver"
    ),
    .before = data
  )

long_track_plot <- long_subtrack_summary |>
  unnest(c(data)) |>
  ggplot() +
  # geom_path(
  #   aes(x = x, y = y, group = particle, color = frame),
  #   linewidth = 0.1
  # ) +
  geom_point(aes(x = x, y = y, color = frame), size = 0.1) +
  facet_wrap(
    nrow = 2,
    facets = vars(date, camera, tissue),
    scales = "free_y"
  ) +
  scale_color_viridis_c() +
  theme_void() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_markdown()
  )
long_track_plot
