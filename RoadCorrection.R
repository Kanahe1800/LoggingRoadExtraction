# load required libs
library(terra)
library(sf)
library(gdistance)
library(Rcpp)
library(raster)

# -------------------------
# User settings (edit here)
# -------------------------
# dir.in  <- "E:/ResRoadsFix/in"
# dir.out <- "E:/ResRoadsFix/out"
# dir.tmp <- "E:/ResRoadsFix/tmp"
## set up directories ##
dir.in <- "C:/Users/mimam/Dropbox/lider_data/input" # for input
dir.out <- "C:/Users/mimam/Dropbox/lider_data/output" # for outputs
dir.tmp <- "C:/Users/mimam/Dropbox/lider_data/tmp" # for temporary files
file_name <- "pine_lake.tif"

dem_file <- file.path(dir.in, file_name)
road_lines <- file.path(dir.in, "road")

# path to your anisotropic C++ source (do NOT change working directory)
anisotropic_cpp <- "C:/Users/sebas/OneDrive/Desktop/Grad school/Scripts/anisotropic_diffusion_2.cpp"


# -------------------------
# Helper functions
# -------------------------

# checks that directory exists
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# checks that file exists
check_file <- function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
}

# normalize to 0-1
rescale01 <- function(x) {
  v <- x[]
  vmin <- min(v, na.rm = TRUE)
  vmax <- max(v, na.rm = TRUE)
  if (is.infinite(vmin) || is.infinite(vmax) || vmax == vmin) {
    x[] <- 0
  } else {
    x[] <- (v - vmin) / (vmax - vmin)
  }
  return(x)
}

# piecewise linear activation 
activation_vec <- function(vals, th, asc = TRUE) {
  stopifnot(length(th) %in% c(1,2))
  start <- if (asc) 0 else 1
  delta <- if (asc) 1 else -1
  
  if (length(th) == 1) {
    min_th <- th[1]
    max_th <- max(vals, na.rm = TRUE)
  } else {
    min_th <- th[1]; max_th <- th[2]
  }
  
  # Avoid division by zero
  if (max_th == min_th) {
    a <- 0
  } else {
    a <- delta / (max_th - min_th)
  }
  b <- start - a * min_th
  y <- a * vals + b
  y[vals > max_th] <- 1 - start
  y[vals < min_th] <- start
  return(round(y, 3))
}

# terra-aware wrapper to apply activation to a raster layer 
activation_raster <- function(r, th, asc = TRUE, filename = "") {
  fun <- function(v) activation_vec(v, th = th, asc = asc)
  app(r, fun, filename = filename)
}

# -------------------------
# Prepare directories & check files
# -------------------------
ensure_dir(dir.out)
ensure_dir(dir.tmp)
check_file(dem_file)
# check_file(points_shp)

if (!file.exists(anisotropic_cpp)) {
  warning("anisotropic C++ file not found at: ", anisotropic_cpp,
          "\nSkipping anisotropic diffusion steps. If you want diffusion, set anisotropic_cpp to the correct path.")
} else {
  sourceCpp(anisotropic_cpp)
}

# -------------------------
# Read DEM and optionally resample
# -------------------------
dem <- rast(dem_file)

# optional: target resampling resolution (meters). Set to NULL to skip resampling.
target_res_m <- NULL

if (!is.null(target_res_m)) {
  # create empty template with same extent / crs and desired resolution
  template <- rast(extent(dem), crs = crs(dem), resolution = target_res_m,
                   nlyrs = 1)
  rs_dem <- resample(dem, template, method = "bilinear")
  dem <- rs_dem
  rm(rs_dem)
}


# -------------------------
# Conductivity parameters
# -------------------------
con_param <- list(
  sigma_min = 0.1,       # to avoid zeros in final surface
  s = c(5, 30),          # slope thresholds
  r = c(0.05, 3),        # roughness thresholds
  g = c(0, 0.6),         # gradient 
  e = 40,                # sobel edge threshold
  alpha = list(r = 1)    # placeholders 
)

# -------------------------
# 1) Slope
# -------------------------
slope <- terrain(dem, v = "slope", unit = "degrees")

sigma_s <- activation_raster(slope, th = con_param$s, asc = FALSE)

# -------------------------
# 2) Roughness
# -------------------------
roughness <- terrain(dem, v = "roughness")

sigma_r <- activation_raster(roughness, th = con_param$r, asc = FALSE)

# -------------------------
# 3) Sobel Edge detection
# -------------------------
sobel_filter_terra <- function(r, size = 3) {
  if (!size %in% c(3,5)) stop("size must be 3 or 5")
  if (size == 3) {
    sobel_x <- matrix(c(-1, 0,  1, -2, 0, 2, -1, 0, 1), nrow=3, byrow=TRUE)
    sobel_y <- matrix(c(-1, -2, -1, 0,0,0, 1,2,1), nrow=3, byrow=TRUE)
  } else {
    sobel_x <- matrix(c(-2, -1,  0,  1,  2,
                        -3, -2,  0,  2,  3,
                        -4, -3,  0,  3,  4,
                        -3, -2,  0,  2,  3,
                        -2, -1,  0,  1,  2), nrow = 5, byrow = TRUE)
    sobel_y <- matrix(c(-2, -3, -4, -3, -2,
                        -1, -2, -3, -2, -1,
                        0,  0,  0,  0,  0,
                        1,  2,  3,  2,  1,
                        2,  3,  4,  3,  2), nrow = 5, byrow = TRUE)
  }
  gx <- focal(r, sobel_x, fun = sum, na.policy = "omit")
  gy <- focal(r, sobel_y, fun = sum, na.policy = "omit")
  sqrt(gx^2 + gy^2)
}

edges <- sobel_filter_terra(slope, size = 5)
sigma_e <- activation_raster(edges, th = con_param$e, asc = FALSE)

# -------------------------
# Combine conductivity surface
# -------------------------
con_surf <- sigma_s * sigma_e * sigma_r

# avoid zeros breaking connectivity
con_surf[is.na(con_surf[])] <- NA
con_surf[con_surf == 0] <- con_param$sigma_min
con_surf <- clamp(con_surf, lower = con_param$sigma_min, upper = Inf)

# -------------------------
# Read roads create endpoints and buffer
# -------------------------
# read road network
roads <- vect(road_lines)

# convert to sf for easier geometry handling
roads_sf <- st_as_sf(roads)

# function to extract start/end points from each line
get_endpoints <- function(line) {
  coords <- st_coordinates(line)
  if (nrow(coords) == 0) return(NULL)
  start_point <- st_point(coords[1, 1:2])
  end_point   <- st_point(coords[nrow(coords), 1:2])
  st_sfc(start_point, end_point, crs = st_crs(line))
}

# apply to each line in the shapefile
endpoints_list <- lapply(st_geometry(roads_sf), get_endpoints)

# combine into sf point layer
endpoints <- st_as_sf(do.call(c, endpoints_list), crs = st_crs(roads_sf))

# assign attributes
endpoints$segment_id <- rep(roads_sf$OBJECTID, each = 2)
endpoints$type <- rep(c("start", "end"), length.out = nrow(endpoints))

# create polygon mask from DEM valid area
mask <- classify(dem, rbind(cbind(NA, NA)), right = NA)  # NA where no data
mask[!is.na(dem)] <- 1
mask_poly <- as.polygons(mask, dissolve = TRUE)

# crop roads by mask extent
roads_vect <- crop(roads, mask_poly)

# buffer roads to restrict least-cost path computation
buff <- buffer(roads_vect, width = 25)
con_surf <- mask(con_surf, buff)

# -------------------------
# Optional anisotropic diffusion / edge enhancement
# Requires anisotropic_diffusion_cpp() available from sourceCpp
# -------------------------
if (exists("anisotropic_diffusion_cpp")) {
  anisotropic_diffusion_filter <- function(x, iteration = 50, lambda = 0.2, k = 10) {
    M <- raster::as.matrix(x)
    M2 <- anisotropic_diffusion_cpp(M * 255, iteration, lambda, k)
    M2 <- M2 / 255
    y <- x
    y[] <- as.vector(M2)
    return(y)
  }
  edge_enhancement <- function(x, iteration = 50, lambda = 0.05, k = 10, pass = 2) {
    for (i in seq_len(pass)) {
      x <- anisotropic_diffusion_filter(x, iteration, lambda, k)
    }
    rescale01(x)
  }
  con_surf <- edge_enhancement(con_surf, iteration = 50, lambda = 0.05, k = 10, pass = 3)
} else {
  message("Anisotropic diffusion skipped: C++ function not found (anisotropic_cpp path).")
}

# -------------------------
# ensure endpoints are ordered and split into start/end pairs
# -------------------------
endpoints <- endpoints[order(endpoints$segment_id, endpoints$type), ]
start_points <- endpoints[endpoints$type == "start", ]
end_points   <- endpoints[endpoints$type == "end", ]

if (nrow(start_points) != nrow(end_points)) {
  warning("Number of start and end points differ. Some pairs may be skipped.")
}

# convert con_surf to raster (raster package) for gdistance transition
con_raster <- raster(con_surf)
rm(con_surf)
gc() # garbage collect

# -------------------------
# Build transition layer and run least-cost paths
# -------------------------
# Use mean transition and correct for diagonal movement
tr <- transition(con_raster, transitionFunction = mean, directions = 8)
tr <- geoCorrection(tr, type = "c")

lcp_list <- vector("list", min(nrow(start_points), nrow(end_points)))
k <- 1
for (i in seq_len(min(nrow(start_points), nrow(end_points)))) {
  s_xy <- st_coordinates(start_points[i, ])
  e_xy <- st_coordinates(end_points[i, ])
  if (any(is.na(s_xy)) || any(is.na(e_xy))) next
  
  # shortestPath can throw an error; catch and continue
  tryCatch({
    sp <- shortestPath(tr, origin = s_xy, goal = e_xy, output = "SpatialLines")
    sp_sf <- st_as_sf(sp)
    sp_sf$segment_id <- start_points$segment_id[i]
    lcp_list[[k]] <- sp_sf
    k <- k + 1
  }, error = function(e) {
    message(sprintf("LCP failed for pair %d -> %d: %s", i, i, e$message))
  })
}

# combine results
lcp_list <- Filter(Negate(is.null), lcp_list)
if (length(lcp_list) == 0) {
  stop("No least-cost paths were computed.")
}
lcp_sf <- do.call(rbind, lcp_list)

# Convert to terra vect for writing if desired
lcp_vect <- vect(lcp_sf)

# -------------------------
# Output
# -------------------------
out_shp <- file.path(dir.out, "RobsonRoad.shp")
writeVector(lcp_vect, out_shp, overwrite = TRUE)