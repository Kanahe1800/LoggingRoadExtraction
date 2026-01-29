library(terra)
library(sf)
library(dbscan)

source("C:/Users/mimam/Dropbox/LoggingRoadExtraction/skeletonize.R")

to_utm_epsg <- function(sfobj_ll) {
  ctd  <- st_coordinates(st_centroid(st_union(sfobj_ll)))
  lon  <- ctd[1]; lat <- ctd[2]
  zone <- floor((lon + 180) / 6) + 1
  if (lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
}

seeds_to_points <- function(
    seeds,
    out_points,
    prob_threshold      = 0.25,
    do_opening          = FALSE,
    opening_size        = 3,
    min_skel_blob_px    = 1,
    min_component_len_m = 1.0,
    min_spacing_m       = 8.0
){
  stopifnot(inherits(seeds, "SpatRaster"))
  
  # Work in meters if input is lon/lat
  if (terra::is.lonlat(seeds)) {
    cent_ll <- st_as_sfc(st_bbox(seeds), crs = 4326) |> st_centroid()
    utm_epsg <- to_utm_epsg(cent_ll)
    message("Reprojecting seeds to ", utm_epsg, " for metric processing")
    seeds <- terra::project(seeds, utm_epsg, method = "bilinear")
  }
  
  # Pixel size (diagonal) in meters, as in np.hypot
  pix_m <- sqrt(sum(terra::res(seeds)^2))
  
  # --- Threshold to binary mask ---
  m <- as.matrix(seeds)
  finite_mat <- is.finite(m)
  if (!any(finite_mat)) stop("No finite pixels in seeds")
  
  thr <- as.numeric(prob_threshold)
  mask_mat <- (m >= thr) & finite_mat
  cat(sprintf("[A] threshold=%.3f → mask pixels: %d\n", thr, sum(mask_mat, na.rm = TRUE)))
  
  # Optional tiny opening (erode→dilate) using focal min/max
  if (do_opening) {
    w <- matrix(1, opening_size, opening_size)
    r_mask <- rast(ext = ext(seeds), crs = crs(seeds), nrow = nrow(seeds), ncol = ncol(seeds))
    values(r_mask) <- as.integer(mask_mat)
    eroded <- terra::focal(r_mask, w = w, fun = min, na.policy = "omit")
    opened <- terra::focal(eroded, w = w, fun = max, na.policy = "omit")
    mask_mat <- as.matrix(opened) > 0.5
    cat(sprintf("[A] after opening: %d pixels\n", sum(mask_mat, na.rm = TRUE)))
  }
  
  # --- Skeletonize (your sourced function must return logical matrix of same dim) ---
  skel_mat <- thin_zhang_suen(mask_mat)
  if (!is.matrix(skel_mat) || any(dim(skel_mat) != dim(mask_mat))) {
    stop("skeletonize.R: thin_zhang_suen() must return a logical/numeric matrix of same dimensions as input.")
  }
  skel_mat <- skel_mat > 0
  cat(sprintf("[B] skeleton pixels: %d\n", sum(skel_mat, na.rm = TRUE)))
  
  # Build skeleton raster with background = NA (CRITICAL for patches())
  skel_r <- rast(ext = ext(seeds), crs = crs(seeds), nrow = nrow(seeds), ncol = ncol(seeds))
  values(skel_r) <- NA_integer_
  skel_r[skel_mat] <- 1L
  
  # --- Remove tiny skeleton fragments (in pixels) ---
  if (min_skel_blob_px > 1) {
    lab <- terra::patches(skel_r, directions = 8)  # labels (NA outside)
    ft  <- as.data.frame(terra::freq(lab))         # <-- no useNA arg
    if (nrow(ft) == 0) {
      cat("[B] no labeled pixels; nothing to keep\n")
      return(invisible(st_sf()))
    }
    keep_ids <- ft$value[ft$count >= min_skel_blob_px]
    skel_keep <- lab %in% keep_ids
  } else {
    skel_keep <- !is.na(skel_r)
  }
  cat(sprintf("[B] skeleton after small-object removal: %d px\n",
              terra::global(skel_keep, "sum", na.rm = TRUE)[[1]]))
  
  # --- Keep long components only (by length in meters) ---
  lab_len <- terra::patches(terra::mask(skel_r, skel_keep, maskvalues = NA), directions = 8)
  props   <- as.data.frame(terra::freq(lab_len))   # <-- no useNA arg
  if (nrow(props) == 0) stop("No skeleton pixels survived. Lower threshold or relax filters.")
  
  props$length_m <- props$count * pix_m
  keep_labels <- props$value[props$length_m >= min_component_len_m]
  
  if (length(keep_labels) == 0) {
    keep_labels <- head(props$value[order(-props$count)], 5)
    message("[C] no component ≥ min length; keeping top 5 longest as fallback.")
  }
  
  # --- Convert logical mask to integer raster with background = NA ---
  skel_final <- (lab_len %in% keep_labels)
  r_final <- terra::ifel(skel_final, 1L, NA_integer_)
  names(r_final) <- "v"
  
  kept_px <- terra::global(r_final, "sum", na.rm = TRUE)[[1]]
  cat(sprintf("[C] kept skeleton pixels (after length filter): %d\n", kept_px))
  if (is.na(kept_px) || kept_px == 0)
    stop("No skeleton pixels survived. Lower threshold or relax filters.")
  
  # --- Raster → points; attach score from seeds raster ---
  xy <- terra::as.data.frame(r_final, xy = TRUE, na.rm = TRUE)
  # xy has columns: x, y, v (v==1 for skeleton pixels)
  
  # sample scores from the original seeds raster at those x,y
  scores <- as.numeric(terra::extract(seeds, xy[, c("x","y")], cells = FALSE)[, 1])
  
  # IMPORTANT: coords expects column NAMES, not the data
  gdf <- sf::st_as_sf(data.frame(score = scores, x = xy$x, y = xy$y),
                      coords = c("x","y"),
                      crs = terra::crs(seeds))
  cat(sprintf("[D] converted to %d skeleton points\n", nrow(gdf)))
  
  # --- Thin points to min_spacing_m (greedy, highest score first) ---
  # choose UTM from geometry centroid
  cent_ll <- sf::st_transform(sf::st_as_sfc(st_bbox(gdf), crs = sf::st_crs(gdf)),
                              4326) |>
    sf::st_centroid()
  utm_epsg <- to_utm_epsg(cent_ll)
  gdf_m <- sf::st_transform(gdf, utm_epsg)
  
  XY <- sf::st_coordinates(gdf_m)
  order_idx <- order(-gdf_m$score)
  alive <- rep(TRUE, nrow(gdf_m))
  keep_idx <- integer(0)
  
  nbr <- dbscan::frNN(XY, eps = min_spacing_m)
  for (i in order_idx) {
    if (!alive[i]) next
    keep_idx <- c(keep_idx, i)
    nbi <- unlist(nbr$id[i])
    if (length(nbi)) alive[nbi] <- FALSE
    alive[i] <- TRUE
  }
  
  gdf_thin <- gdf[keep_idx, ]
  cat(sprintf("[E] thinned points: %d\n", nrow(gdf_thin)))
  gdf_thin <- sf::st_transform(gdf_thin, 4326)
  
  # --- Save ---
  sf::st_write(gdf_thin, out_points, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  cat(sprintf("Saved %d points → %s\n", nrow(gdf_thin), out_points))
  
}


out_points <- "C:/Users/mimam/Dropbox/lider_data/output/seed_points_test_v2.geojson"

pts <- seeds_to_points(
  seeds = seeds,
  out_points = out_points,
  prob_threshold = 0.4,
  min_spacing_m = 10.0,           # Adjust this to control point density
  min_component_len_m = 2.0
)