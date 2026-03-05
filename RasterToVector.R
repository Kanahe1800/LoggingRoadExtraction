library(terra)
library(sf)
library(dbscan)
library(RANN)

source("EnvLoader.R")
source(skeletonize_path)

to_utm_epsg <- function(sfobj_ll) {
  # 1. Get centroid coordinates
  ctd  <- sf::st_coordinates(sf::st_centroid(sf::st_union(sfobj_ll)))
  lon  <- ctd[1]; lat <- ctd[2]
  
  # 2. Calculate UTM Zone (1-60)
  zone <- floor((lon + 180) / 6) + 1
  
  # 3. Validate Zone (Must be 1-60)
  if (zone < 1 || zone > 60) {
    warning("Calculated UTM zone is invalid. Defaulting to WGS84 (EPSG:4326)")
    return("EPSG:4326")
  }
  
  # 4. Construct EPSG code
  epsg <- if (lat >= 0) (32600 + zone) else (32700 + zone)
  return(paste0("EPSG:", epsg))
}

seeds_to_points <- function(
    seeds,
    out_points,
    prob_threshold,
    do_opening,
    opening_size,
    min_skel_blob_px,
    min_component_len_m,
    min_spacing_m,
    snap_radius  # Search radius (in pixels) to find the true road center
){
  stopifnot(inherits(seeds, "SpatRaster"))
  
  # 1. Coordinate System Check
  if (terra::is.lonlat(seeds)) {
    cent_ll <- st_as_sfc(st_bbox(seeds), crs = 4326) |> st_centroid()
    utm_epsg <- to_utm_epsg(cent_ll)
    seeds <- terra::project(seeds, utm_epsg, method = "bilinear")
  }
  
  pix_m <- sqrt(sum(terra::res(seeds)^2))
  
  # 2. Thresholding and Skeletonization
  m <- as.matrix(seeds)
  mask_mat <- (m >= as.numeric(prob_threshold)) & is.finite(m)
  
  if (do_opening) {
    w <- matrix(1, opening_size, opening_size)
    r_mask <- rast(seeds); values(r_mask) <- as.integer(mask_mat)
    mask_mat <- as.matrix(terra::focal(terra::focal(r_mask, w, min), w, max)) > 0.5
  }
  
  skel_mat <- thin_zhang_suen(mask_mat) > 0
  skel_r   <- rast(seeds); values(skel_r) <- NA_integer_; skel_r[skel_mat] <- 1L
  
  # 3. Component Filtering (Size and Length)
  lab <- terra::patches(skel_r, directions = 8)
  if (min_skel_blob_px > 1) {
    ft <- as.data.frame(terra::freq(lab))
    skel_r <- terra::ifel(lab %in% ft$value[ft$count >= min_skel_blob_px], 1L, NA_integer_)
  }
  
  lab_len <- terra::patches(skel_r, directions = 8)
  props   <- as.data.frame(terra::freq(lab_len))
  props$length_m <- props$count * pix_m
  r_final <- terra::ifel(lab_len %in% props$value[props$length_m >= min_component_len_m], 1L, NA_integer_)
  
  # 4. FAST GLOBAL SNAPPING (Focal Method)
  message(sprintf("Processing points with Global Focal Snapping..."))
  
  # A. Extract skeleton coordinates FIRST (before snapping)
  xy_skel <- terra::as.data.frame(r_final, xy = TRUE, na.rm = TRUE)
  
  # B. Calculate local maxima for the entire raster
  size <- (snap_radius * 2) + 1
  focal_max <- terra::focal(seeds, w = size, fun = "max", na.rm = TRUE)
  is_local_max <- (seeds == focal_max)
  
  # C. Mask maxima to skeleton vicinity to keep RAM usage low
  search_dist <- (snap_radius + 1) * pix_m
  skel_poly <- sf::st_as_sf(terra::as.polygons(r_final))
  skel_buffer <- sf::st_buffer(skel_poly, dist = search_dist)
  
  local_max_cropped <- terra::mask(is_local_max, terra::vect(skel_buffer))
  max_pts <- terra::as.data.frame(local_max_cropped, xy = TRUE, na.rm = TRUE)
  max_pts <- max_pts[max_pts[[3]] == 1, 1:2] 
  
  # D. Snap: Find nearest high-probability peak for each skeleton point
  if (nrow(max_pts) > 0) {
    tree <- nn2(data = max_pts, query = xy_skel[,1:2], k = 1)
    xy_snapped <- max_pts[tree$nn.idx, ]
  } else {
    xy_snapped <- xy_skel[,1:2] # Fallback if no peaks found
  }
  
  # 5. CONVERT TO SF AND THIN
  message("Creating spatial features and thinning points...")
  
  # Ensure we have snapped coordinates
  if (!exists("xy_snapped") || nrow(xy_snapped) == 0) {
    stop("Snapping failed: No valid road peaks found. Try lowering 'prob_threshold'.")
  }
  
  # A. Create the 'gdf' object FIRST
  scores <- as.numeric(terra::extract(seeds, xy_snapped)[,1])
  gdf <- sf::st_as_sf(data.frame(score = scores, x = xy_snapped[,1], y = xy_snapped[,2]),
                      coords = c("x","y"), crs = terra::crs(seeds))
  
  # B. Prepare for thinning (Project to UTM for accurate distance in meters)
  # Use a fixed, correct EPSG for your region (Oak Bay, BC = 32610)
  target_epsg <- 32610 
  gdf_m <- sf::st_transform(gdf, target_epsg)
  
  # C. GREEDY THINNING (Grid-Based Optimization)
  print("Final Greedy Thinning using Grid...")
  
  # 1. Map points to a grid based on min_spacing_m
  # Using gdf_m (the projected version) ensures units are in meters
  print("1")
  coords <- sf::st_coordinates(gdf_m)
  grid_x <- floor(coords[,1] / min_spacing_m)
  grid_y <- floor(coords[,2] / min_spacing_m)
  
  # 2. Assign each point to a cell ID
  print("2")
  cell_id <- paste(grid_x, grid_y, sep = "_")
  
  # 3. Create a data frame to identify the best point per cell
  print("3")
  df_thin <- data.frame(
    id = 1:nrow(gdf_m),
    cell_id = cell_id,
    score = gdf_m$score
  )
  
  # 4. Pick the index of the point with the max score per cell
  # We use aggregate to find the row index in gdf_m (and gdf)
  print("4")
  best_indices <- aggregate(id ~ cell_id, data = df_thin, FUN = function(x) {
    x[which.max(df_thin$score[x])]
  })$id
  
  # 5. Subset the original (unprojected) gdf using those indices
  print("5")
  gdf_thin <- gdf[best_indices, ]
  
  # 6. Final Write
  print("6")
  sf::st_write(gdf_thin, out_points, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  
  return(gdf_thin)
}


pts <- seeds_to_points(
  seeds = seeds,
  out_points_file,
  prob_threshold      = prob_threshold,
  do_opening          = do_opening,
  opening_size        = opening_size,
  min_skel_blob_px    = min_skel_blob_px,
  min_component_len_m = min_component_len_m,
  min_spacing_m       = min_spacing_m,
  snap_radius = 4
)