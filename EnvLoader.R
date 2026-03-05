# ══════════════════════════════════════════════════════════════════════════════
# env_loader.R — source this at the top of every script
# Usage: source("env_loader.R")
# ══════════════════════════════════════════════════════════════════════════════

# ── Load both env files ────────────────────────────────────────────────────────
if (!file.exists("local.env")) {
  stop("local.env not found. Copy local.env.example → local.env and fill in your paths.")
}
readRenviron("local.env")
readRenviron("settings.env")

# ── Paths (from local.env) ─────────────────────────────────────────────────────
dir.in           <- Sys.getenv("DIR_IN")
dir.out          <- Sys.getenv("DIR_OUT")
dir.tmp          <- Sys.getenv("DIR_TMP")
skeletonize_path <- Sys.getenv("SKELETONIZE_PATH")

# ── File names (from settings.env) ────────────────────────────────────────────
file_name             <- Sys.getenv("FILE_NAME")
seeds_file            <- file.path(dir.out, Sys.getenv("SEEDS_FILE"))
out_points_file       <- file.path(dir.out, Sys.getenv("OUT_POINTS_FILE"))
edited_points_file    <- file.path(dir.out, Sys.getenv("EDITED_FILE"))
final_road_path_file  <- file.path(dir.out, Sys.getenv("FINAL_ROAD_LINE"))


# ── createSeeds parameters ────────────────────────────────────────────────────
slope_smooth_window <- as.integer(Sys.getenv("SLOPE_SMOOTH_WINDOW"))
roughness_sd_small  <- as.integer(Sys.getenv("ROUGHNESS_SD_SMALL"))
roughness_sd_med    <- as.integer(Sys.getenv("ROUGHNESS_SD_MED"))
roughness_sd_large  <- as.integer(Sys.getenv("ROUGHNESS_SD_LARGE"))
sobel_size          <- as.integer(Sys.getenv("SOBEL_SIZE"))
gaussian_size       <- as.integer(Sys.getenv("GAUSSIAN_SIZE"))
gaussian_sigma      <- as.numeric(Sys.getenv("GAUSSIAN_SIGMA"))

# ── RasterToVector parameters ─────────────────────────────────────────────────
prob_threshold      <- as.numeric(Sys.getenv("PROB_THRESHOLD"))
do_opening          <- as.logical(Sys.getenv("DO_OPENING"))
opening_size        <- as.integer(Sys.getenv("OPENING_SIZE"))
min_skel_blob_px    <- as.integer(Sys.getenv("MIN_SKEL_BLOB_PX"))
min_component_len_m <- as.numeric(Sys.getenv("MIN_COMPONENT_LEN_M"))
min_spacing_m       <- as.numeric(Sys.getenv("MIN_SPACING_M"))

# ── RoadExtraction parameters ─────────────────────────────────────────────────
cost_barrier_threshold_small  <- as.numeric(Sys.getenv("COST_BARRIER_THRESHOLD_SMALL"))
cost_barrier_threshold_big    <- as.numeric(Sys.getenv("COST_BARRIER_THRESHOLD_BIG"))
cost_barrier_value            <- as.numeric(Sys.getenv("COST_BARRIER_VALUE"))
tr_directions                 <- as.integer(Sys.getenv("TR_DIRECTIONS"))
max_distance_threshold        <- as.numeric(Sys.getenv("MAX_DISTANCE_THRESHOLD"))
mean_cost_threshold           <- as.numeric(Sys.getenv("MEAN_COST_THRESHOLD"))
min_segments_per_line         <- as.integer(Sys.getenv("MIN_SEGMENTS_PER_LINE"))

# ── Shiny parameters ──────────────────────────────────────────────────────────
shiny_default_lng  <- as.numeric(Sys.getenv("SHINY_DEFAULT_LNG"))
shiny_default_lat  <- as.numeric(Sys.getenv("SHINY_DEFAULT_LAT"))
shiny_default_zoom <- as.integer(Sys.getenv("SHINY_DEFAULT_ZOOM"))

message("✓ Environment loaded from local.env + settings.env")