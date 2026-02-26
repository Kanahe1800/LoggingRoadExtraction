### road extraction 2.0 ###

## required libs ##
library(terra)
library(whitebox)
library(autothresholdr)

## set up directories ##
dir.in <- "E:/MSc/LiDAR/ChipmunkCreek/" # for inputs
dir.out <- "E:/ResRoadsFix/out/" # for outputs
dir.tmp <- "E:/ResRoadsFix/tmp/" # for temporary files

## custom functions ##
# activation function
activation <- function(x, th, asc = TRUE) {
  stopifnot(length(th) %in% c(1, 2))  # Ensure 1 or 2 threshold values
  
  # Extract raster values
  x <- x[]
  
  # Define if ascending or descending
  start <- if (asc) 0 else 1
  delta <- if (asc) 1 else -1
  
  if (length(th) == 1) {
    # If only min threshold is provided, scale from `th[1]` to max(x)
    min_th <- th[1]
    max_th <- max(x, na.rm = TRUE)  # Use max value in raster
    
  } else {
    # If two values are provided, use them as min and max thresholds
    min_th <- th[1]
    max_th <- th[2]
  }
  
  # Compute slope (a) and intercept (b) for the linear equation
  a <- delta / (max_th - min_th)
  b <- start - a * min_th
  
  # Apply linear transformation
  y <- a * x + b
  
  # Restrict values to the threshold range
  y[x > max_th] <- 1 - start
  y[x < min_th] <- start
  
  return(round(y, 3))
}

# sobel edge detection
sobel_filter <- function(r, size = 3, output = c("magnitude", "direction", "both")) {
  output <- match.arg(output)  # ensures valid choice
  
  if (!size %in% c(3, 5)) {
    stop("Size must be either 3 or 5.")
  }
  
  # Define Sobel kernels
  if (size == 3) {
    sobel_x <- matrix(c(-1, 0,  1,
                        -2, 0,  2,
                        -1, 0,  1), nrow=3, byrow=TRUE)
    sobel_y <- matrix(c(-1, -2, -1,
                        0,  0,  0,
                        1,  2,  1), nrow=3, byrow=TRUE)
  } else {
    sobel_x <- matrix(c(-2, -1,  0,  1,  2,
                        -3, -2,  0,  2,  3,
                        -4, -3,  0,  3,  4,
                        -3, -2,  0,  2,  3,
                        -2, -1,  0,  1,  2), nrow=5, byrow=TRUE)
    sobel_y <- matrix(c(-2, -3, -4, -3, -2,
                        -1, -2, -3, -2, -1,
                        0,  0,  0,  0,  0,
                        1,  2,  3,  2,  1,
                        2,  3,  4,  3,  2), nrow=5, byrow=TRUE)
  }
  
  # Apply focal filters
  gx <- focal(r, w = sobel_x, fun = sum, na.policy = "omit")
  gy <- focal(r, w = sobel_y, fun = sum, na.policy = "omit")
  
  # Compute magnitude and direction
  mag <- sqrt(gx^2 + gy^2)
  dir <- atan2(gy, gx)
  
  # Return based on user choice
  if (output == "magnitude") return(mag)
  if (output == "direction") return(dir)
  if (output == "both") return(list(magnitude = mag, direction = dir))
}

# gaussian filter
gaussian_kernel <- function(size = 5, sigma = 1) {
  center <- floor(size / 2)
  x <- -center:center
  y <- -center:center
  kernel <- outer(x, y, function(x, y) exp(-(x^2 + y^2) / (2 * sigma^2)))
  kernel / sum(kernel)  # Normalize so the sum is 1
}

# get Otsu threshold for raster
otsu_thresh <- function(r) {
  vals <- values(r, na.rm = TRUE)
  min <- min(vals)
  max <- max(vals)
  
  vals_int <- as.integer(255 * (vals - min) / (max - min))
  
  thr <- auto_thresh(vals_int, method = "Otsu")
  
  thr_real <- min + (thr / 255) * (max - min)
  
  return(thr_real[1])
}

#### terrain variables ####
## call in the DEM ##
dem <-  terra::rast(file.path(dir.in, "ChipmunkCreek_DEM_Clip.tif"))


## slope (difference) ##
# get slope
slope <- terrain(dem, "slope")

# make a low pass filter window
w <- matrix(1, 51, 51) # uniform weighting
w <- w / sum(w)        # normalize

# apply low pass filter to the slope layer
slope_smooth <- focal(slope, w = w, fun = "sum", na.policy = "omit")

# subtract the smoothed slope from the raw slope
slope_dif <- slope - slope_smooth

# clean memory
rm(slope_smooth)
rm(w)
gc()


## multiscale roughness ##
# run sd filter on slope layer at different window sizes (3, 7, 15 or 5, 11, 21)
sd5  <- focal(slope, w = matrix(1,5,5), fun = sd, na.rm=TRUE)
sd11  <- focal(slope, w = matrix(1,11,11), fun = sd, na.rm=TRUE)
sd21 <- focal(slope, w = matrix(1,21,21), fun = sd, na.rm=TRUE)

# combine scales: roads should have low roughness at multiple scales
roughness_ms <- (sd5 + sd11 + sd21) / 3

# clean memory
rm(sd5)
rm(sd11)
rm(sd21)
gc()


## profile curvature ##
# calculate profile curvature layer
wbt_profile_curvature(
  dem = file.path(dir.in, "ChipmunkCreek_DEM_Clip.tif"),
  output = file.path(dir.out, "profcurv_cc.tif")
)

# call the profile curvature layer in
prof_curv <- rast(file.path(dir.out, "profcurv_cc.tif"))

# apply sobel edge detection
edges <- sobel_filter(prof_curv, size = 5, output = "magnitude")
edges <- ifel(is.na(edges), 0.01, edges) # remove some introduced NA's

# build a gaussian filter kernel
kernel <- gaussian_kernel(size = 21, sigma = 3)

# apply gaussian filter
edges <- focal(edges, w = kernel, fun = sum, na.policy = "omit", pad = TRUE)
edges <- ifel(edges < 0.01, NA, edges) # remove the NA padding from previous step

# clean memory
rm(prof_curv)
rm(kernel)
gc()

## Otsu filtering ##

# multiscale roughness
thr_rough <- otsu_thresh(roughness_ms)
plot(roughness_ms > thr_rough, main = "roughness")

# slope (dif) # this doesnt work so well
thr_slope <- otsu_thresh(slope_dif)
plot(slope_dif < thr_slope, main = "slope")

# edges
thr_edge <- otsu_thresh(edges)
plot(edges > thr_edge)

# use overlap of Otsu thresholds to get seed pixels
seeds <- ifel(edges > thr_edge, 0.333, 0) + ifel(roughness_ms > thr_rough, 0.333, 0) + ifel(slope_dif < thr_slope, 0.333, 0)
plot(seeds)

## Converting seed layer to seed points ##

# use a threshold to define seed pixels
seed_candidates <- ifel(seeds > 0.9, 1, NA)

###############################################################################
# potentially resample to a coarser resolution to make fewer candidate points #
empty_rast <- rast(ext = ext(seeds), crs = crs(seeds), res = 100) # eg. resolution = 100m
seed_candidates <- resample(seeds, empty_rast)
###############################################################################

# convert the seed candidate pixels to points
seed_pts <- as.points(seed_candidates, values = FALSE, na.rm = TRUE)

####################################
# for use with 'shiny' application #
library(sf)

# convert from spatvector object to sf object
seed_sf <- st_as_sf(seeds_pts)

# Reproject to WGS84 lon/lat for using a base map in the application
seed_sf <- st_transform(seed_sf, crs = 4326)
####################################