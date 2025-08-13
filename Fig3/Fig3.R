library(tidyverse)
library(ggtext)
library(ggbeeswarm)
library(cowplot)
library(here)

source(here("utils", "helper_functions.R"))


# import ----------------------------------------------------------------
# start at line 27

# df <- feathers |>
#   # only keep first hour of tracking
#   filter(frame < 3600 * 15)

# # keep any track > 5 seconds (15 FPS)
# frame_filter <- df |>
#   group_by(video, particle) |>
#   tally() |>
#   filter(n > 75)

# filtered <- frame_filter |>
#   select(-n) |>
#   left_join(df)

filtered <- read_rds(here("Fig3", "filtered.rds"))

# nest_cols <- c("date", "video", "particle")
# nested <- quick_nest(filtered, nest_cols)

# track_summary <- calculate_track_features_parallel(
#   nested,
#   fps = 15,
#   chunk_size = 500
# )

track_summary <- read_rds(here("Fig3", "track_summary.rds"))

half_arena_tracks <- track_summary |>
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
    net_displacement > 50,
    sd_x > 5 & sd_y > 5
  )

# plot the tracks
track_plot <- half_arena_tracks |>
  # remove vidoes without tracks
  filter(case_when(
    date == '20250529' & camera == '24568744' ~ FALSE,
    date == '20250626' & camera == '25112214' ~ FALSE,
    date == '20250703' & camera == '24568744' ~ FALSE,
    date == '20250703' & camera == '25112214' ~ FALSE,
    date == '20250717' & camera == '24568709' ~ FALSE,
    TRUE ~ TRUE
  )) |>
  ggplot() +
  geom_path(
    aes(x = x, y = y, group = particle, color = frame),
    linewidth = 0.1
  ) +
  facet_wrap(
    nrow = 2,
    facets = vars(date, camera, tissue),
    scales = "free_y"
  ) +
  scale_color_viridis_c(option = 'mako') +
  scale_x_continuous(
    breaks = seq(
      min(half_arena_tracks$x),
      max(half_arena_tracks$x),
      length.out = 5
    ),
    expand = c(0, 0),
    labels = c("0", "9.125", "18.25", "27.375", "36.5")
  ) +
  theme_void() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_markdown()
  )

save_plot(
  here("Fig3", "SuppFig1.pdf"),
  track_plot,
  base_height = 8,
  base_width = 12
)

# subtrack analysis ------------------------------------------------------

# chunked_data <- filtered |>
#   group_by(video, particle) |>
#   group_split() |>
#   map_dfr(~ split_trajectory(.x, frame_rate = 15, chunk_duration_sec = 5))

# nest_cols <- c("date", "video", "particle", "subparticle")
# nested <- quick_nest(chunked_data, nest_cols)

# subtrack_summary <- calculate_track_features_parallel(
#   nested,
#   fps = 15,
#   chunk_size = 1000
# )

subtrack_summary <- read_rds(here("Fig3", 'subtrack_summary.rds')) |>
  mutate(
    date = str_extract(video, "2025[0-9]{4}"),
    camera = str_extract(video, "[0-9]{8}$"),
    # response (generic variable so we can use fit_model across experiments) == region of interest
    tissue = case_when(
      camera %in% c("24568709", "24568744") ~ "intestine",
      camera %in% c("25128038", "25112214") ~ "liver"
    ),
    .before = data
  ) |>
  filter(case_when(
    date == '20250529' & camera == '24568744' ~ FALSE,
    date == '20250626' & camera == '25112214' ~ FALSE,
    date == '20250703' & camera == '24568744' ~ FALSE,
    date == '20250703' & camera == '25112214' ~ FALSE,
    date == '20250717' & camera == '24568709' ~ FALSE,
    TRUE ~ TRUE
  ))

source("~/GitHub/invision-tools/R-functions/model_utils.R")

feature_cols <- subtrack_summary %>%
  ungroup() |>
  select(-(video:sd_y)) %>%
  select_if(is.numeric) %>%
  names()

# model_results_list <- map(
#   feature_cols,
#   ~ fit_model_interaction(
#     .x,
#     subtrack_summary,
#     fixed_effects = c("tissue", "frame_start"),
#     use_temporal_correlation = TRUE
#   )
# )

# names(model_results_list) <- feature_cols

# # Extract the results
# results <- map_dfr(model_results_list, ~ .x$results) |>
#   filter(str_detect(term, "tissue"))

results <- read_rds(here("Fig3", "model_results.rds"))

feature_summary <- subtrack_summary |>
  pivot_longer(
    cols = frame_start:curv_q90,
    names_to = "feature",
    values_to = "value"
  ) |>
  group_by(feature) |>
  filter(is.finite(value)) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

results_prep <- results |>
  filter(!term %in% c('frame_start', "(Intercept)")) |>
  mutate(p_adj = p.adjust(p.value, method = "fdr")) |>
  mutate(
    sig = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ''
    )
  ) |>
  left_join(feature_summary |> select(feature, sd)) |>
  mutate(
    estimate = case_when(
      # interaction effect scaled to per 10 minutes
      str_detect(term, ":") ~ estimate * 15 * 60 * 20,
      TRUE ~ estimate
    ),
    cohens_d = estimate / sd,
    se_cohens_d = std.error / sd
  ) |>
  arrange(desc(cohens_d))

feature_order <- results_prep |>
  group_by(feature) |>
  summarise(mean_cohens_d = mean(cohens_d, na.rm = TRUE), .groups = 'drop') |>
  arrange(mean_cohens_d) |>
  pull(feature)

feature_labels <- c(
  "directional_persistence" = "Directional persistence",
  "heading_autocorr" = "Heading autocorrelation",
  "angvel_q10" = "Angular velocity (10th percentile)",
  "accel_q10" = "Acceleration (10th percentile)",
  "angacc_q10" = "Angular acceleration (10th percentile)",
  "convex_hull_area" = "Convex hull area",
  "speed_autocorr" = "Speed autocorrelation",
  "net_displacement" = "Net displacement",
  "angvel_q50" = "Angular velocity (50th percentile)",
  "speed_q10" = "Speed (10th percentile)",
  "angacc_q50" = "Angular acceleration (50th percentile)",
  "movement_efficiency" = "Movement efficiency",
  "angular_acceleration_mean" = "Mean angular acceleration",
  "jerk_mean" = "Mean jerk",
  "radius_of_gyration" = "Radius of gyration",
  "acceleration_mean" = "Mean acceleration",
  "bbox_aspect_ratio" = "Bounding box aspect ratio",
  "angular_velocity_mean" = "Mean angular velocity",
  "turning_bias" = "Turning bias",
  "accel_q50" = "Acceleration (50th percentile)",
  "mean_heading" = "Mean heading",
  "speed_mean" = "Mean speed",
  "speed_q50" = "Speed (50th percentile)",
  "curv_q10" = "Curvature (10th percentile)",
  "path_length" = "Path length",
  "curv_q50" = "Curvature (50th percentile)",
  "fractal_dimension" = "Fractal dimension",
  "mean_curvature" = "Mean curvature",
  "curv_q90" = "Curvature (90th percentile)",
  "speed_q90" = "Speed (90th percentile)",
  "speed_var" = "Speed variance",
  "acceleration_var" = "Acceleration variance",
  "sinuosity" = "Sinuosity",
  "tortuosity" = "Tortuosity",
  "angular_velocity_var" = "Angular velocity variance",
  "angacc_q90" = "Angular acceleration (90th percentile)",
  "total_turn" = "Total turn",
  "angvel_q90" = "Angular velocity (90th percentile)",
  "accel_q90" = "Acceleration (90th percentile)",
  "heading_variance" = "Heading variance",
  "max_dist_from_start" = "Max distance from start",
  "speed_max" = "Max speed",
  "straightness" = "Straightness"
)

(results_plot <- results_prep |>
  mutate(
    feature = factor(feature, levels = feature_order),
    ,
    term = case_when(
      str_detect(term, ":") ~ "Tissue * Frame",
      TRUE ~ "Tissue"
    ),
  ) |>
  ggplot() +
  geom_tile(
    aes(x = term, y = feature, fill = cohens_d)
  ) +
  geom_text(
    aes(x = term, y = feature, label = sig),
    color = "black",
    size = 2.5,
    vjust = 0.75,
  ) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = feature_labels) +
  scale_fill_distiller(
    limits = c(-1, 1),
    palette = "BrBG",
    breaks = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 5,
      barheight = 1
    )
  ) +
  labs(
    x = "Fixed effect",
    y = "Feature",
    fill = "Effect size"
  ) +
  theme_half_open() +
  theme(
    axis.title = element_markdown(size = 9),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 7),
    legend.title = element_text(size = 9, hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.title.position = "top",
    # plot.margin = margin(5, 5, 5, 20),
    # legend.box.margin = margin(0, 0, 0, -120),
    strip.text = element_markdown(size = 8),
  ) +
  NULL)

(feature_plot <- subtrack_summary |>
  ungroup() |>
  select(
    date,
    tissue,
    frame_start,
    speed_mean,
    heading_variance
  ) |>
  pivot_longer(
    -c(date, tissue, frame_start),
    names_to = "feature",
    values_to = "value"
  ) |>
  mutate(
    feature = case_when(
      feature == "speed_mean" ~ "Mean speed", # fixed diff, liver higher
      feature == "heading_variance" ~ "Heading variance", # fixed diff, intestine higher
      TRUE ~ feature
    ),
    value = case_when(
      feature == "speed_mean" ~ value / 126.5, # scale to mm/s
      TRUE ~ value
    )
  ) |>
  ggplot(aes(x = frame_start / 15 / 60, y = value, color = tissue)) +
  # geom_point(size = 0.1, alpha = 0.1) +
  geom_smooth() +
  facet_grid(rows = vars(feature), scales = "free_y") +
  labs(x = 'Minute', y = 'Value', color = "Tissue source") +
  scale_color_manual(
    values = c("intestine" = "indianred", "liver" = "steelblue"),
    labels = c("intestine" = "Intestine", "liver" = "Liver")
  ) +
  theme_half_open() +
  theme(
    axis.title = element_markdown(size = 9),
    axis.text = element_markdown(size = 8),
    legend.title = element_text(size = 9, hjust = 0.5),
    legend.text = element_text(size = 8),
    # legend.position = "bottom",
    # legend.title.position = "top",
    strip.text = element_markdown(size = 8),
  ) +
  NULL)

## behavioral profiles ------------------------------------------------------------

sample_summary <- subtrack_summary %>%
  group_by(date, video) %>%
  summarise(
    n_tracks = n_distinct(particle),
    tissue = first(tissue),
    .groups = "drop"
  ) %>%
  mutate(sample_id = paste(date, video, sep = "_"))

create_biological_sample_profiles <- function(data, feature_cols) {
  cat("Creating biological sample profiles (date-tissue combinations)...\n")

  # Simple aggregation across technical replicates (videos)
  # This preserves the biological sample-to-sample variation we want for PCA
  sample_profiles <- data %>%
    group_by(date, tissue) %>% # Biological samples = date-tissue combinations
    summarise(
      across(
        all_of(feature_cols),
        ~ {
          # Filter out -Inf, Inf, and NA values before calculating mean
          clean_values <- .x[is.finite(.x) & !is.na(.x)]
          if (length(clean_values) > 0) {
            mean(clean_values)
          } else {
            NA_real_
          }
        }
      ),
      n_videos = n_distinct(video),
      n_tracks = n_distinct(particle),
      .groups = "drop"
    ) %>%
    mutate(sample_id = paste(date, tissue, sep = "_"))

  cat("Biological sample summary:\n")
  print(sample_profiles %>% select(date, tissue, sample_id, n_videos, n_tracks))

  # Check for variation within tissue types
  variation_check <- sample_profiles %>%
    group_by(tissue) %>%
    summarise(
      across(all_of(feature_cols[1:3]), ~ sd(.x, na.rm = TRUE)), # Check first 3 features
      .groups = "drop"
    )

  cat("\nVariation within tissue types (first 3 features):\n")
  print(variation_check)

  return(sample_profiles)
}

sample_profiles <- create_biological_sample_profiles(subtrack_summary, feature_cols)

behavioral_matrix <- sample_profiles %>%
  select(all_of(feature_cols)) %>%
  as.matrix()

rownames(behavioral_matrix) <- sample_profiles$sample_id

pca_result <- prcomp(behavioral_matrix, center = TRUE, scale. = TRUE)

pca_coords <- tibble(
  sample_id = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = if (ncol(pca_result$x) >= 3) pca_result$x[, 3] else NA
) %>%
  separate(sample_id, into = c("date", "tissue"), sep = "_")

var_explained <- summary(pca_result)$importance[2, ] * 100

(pca_plot <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(aes(shape = date, fill = tissue), size = 4, show.legend = FALSE) +
  # geom_text(aes(label = date)) +
  stat_ellipse(type = "norm", level = 0.68) + # 1 SD
  scale_color_manual(
    values = c("intestine" = "indianred", "liver" = "steelblue"),
    labels = c("intestine" = "Intestine", "liver" = "Liver")
  ) +
  scale_fill_manual(
    values = c("intestine" = "indianred", "liver" = "steelblue"),
    labels = c("intestine" = "Intestine", "liver" = "Liver")
  ) +
  scale_shape_manual(values = c(15, 16, 17, 18, 25)) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
    color = "Tissue source",
    shape = "Replicate"
  ) +
  theme_half_open() +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_markdown(size = 8),
    legend.position = 'empty'
  ) +
  NULL)


example_tracks <- filtered |>
  filter(
    date == "20250717",
    str_detect(video, "24568744|25128038")
  ) |>
  mutate(
    camera = str_extract(video, "[0-9]{8}$"),
    tissue = case_when(
      camera %in% c("24568709", "24568744") ~ "intestine",
      camera %in% c("25128038", "25112214") ~ "liver"
    ),
    .before = particle
  )

example_summary <- example_tracks |>
  group_by(camera, tissue, particle) |>
  arrange(frame) |>
  summarise(
    n_frames = n(),
    chord = sqrt((last(x) - first(x))^2 + (last(y) - first(y))^2)
  ) |>
  group_by(camera, tissue) |>
  slice_max(chord, n = 75)

(ex_track_plot <- example_summary |>
  left_join(example_tracks) |>
  left_join(
    distinct(select(
      ungroup(subtrack_summary),
      video,
      particle,
      tissue,
      frame = frame_start,
      speed_mean,
      heading_variance
    ))
  ) |>
  fill(speed_mean, heading_variance, .direction = "down") |>
  ggplot() +
  geom_path(
    aes(x = x, y = y, group = particle, color = speed_mean / 126.5), # 126.5 px/s
    linewidth = 1
  ) +
  scico::scale_color_scico(
    palette = "batlow",
  ) +
  scale_x_continuous(
    breaks = seq(min(filtered$x), max(filtered$x), length.out = 5),
    expand = c(0, 0),
    labels = c("", "6.75", "13.5", "20.25", "27")
  ) +
  scale_y_continuous(
    breaks = seq(min(filtered$y), max(filtered$y), length.out = 5),
    expand = c(0, 0),
    labels = c("", "6.75", "13.5", "20.25", "27")
  ) +
  facet_wrap(
    nrow = 2,
    facets = vars(tissue),
    scales = "free",
    labeller = labeller(tissue = c("intestine" = "Intestine", "liver" = "Liver"))
  ) +
  labs(x = "X (mm)", y = "Y (mm)", color = "Mean speed (mm/s)") +
  theme_half_open() +
  theme(
    aspect.ratio = 1,
    axis.title = element_markdown(size = 9),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 7.5),
    legend.title = element_text(size = 9, hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.title.position = "top",
    # plot.margin = margin(5, 5, 5, 20),
    # legend.box.margin = margin(0, 0, 0, -120),
    strip.text = element_markdown(size = 8),
  ) +
  NULL)

top <- plot_grid(
  results_plot,
  ex_track_plot,
  ncol = 2,
  align = "h",
  axis = 'b',
  rel_widths = c(1.4, 1),
  labels = c("A", "B")
)

bottom <- plot_grid(
  pca_plot,
  feature_plot,
  ncol = 2,
  align = "h",
  axis = "b",
  rel_widths = c(1, 2),
  labels = c("C", "D")
)

final_plot <- plot_grid(
  top,
  bottom,
  ncol = 1,
  rel_heights = c(2, 1.25)
)

save_plot(
  here("Fig3", "Fig3.pdf"),
  final_plot,
  base_width = 6.5,
  base_height = 8
)

save_plot(
  here("Fig3", "Fig3.png"),
  final_plot,
  base_width = 6.5,
  base_height = 8,
  bg = 'white'
)
