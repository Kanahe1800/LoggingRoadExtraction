## SURFACE DETRENDING TRANSLATED TO PYTHON ##
## translated for MATLAB script: https://github.com/HydrogeomorphologyTools/DTM-Inpainting-surface-roughness-restitution/tree/master
## which applies methods used in:  https://doi-org.ezproxy.library.uvic.ca/10.1002/esp.4739

################
### COMMENTS ###
################

#* Currently this function is very slow, need to optimize
#* need to look into if it can be sped up by restricting the area used to define
#* the residuals i.e., variance, e.g., restrict to 82m buffer?
#* also need to look in to if using non float 64 can speed up processing.
#* 
#* in the R wrapper I should add a bit that converts an input raster to the proper
#* matrix type and then back to a raster using the original raster ext and crs

############################
### INSTALL DEPENDANCIES ###
############################

# install.packages("reticulate") # if required
# install.packages("terra") # if required

library(reticulate)
library(terra)

###########################
### DEFINE THE FUNCTION ###
###########################

py_run_string("
import numpy as np
from scipy.ndimage import convolve, generic_filter
import random

def detrend_surface(Z, ker_size=9, big_window=41, n_samples=1):
    # Convert Z and A to numpy arrays (if they are not already)
    Z = np.array(Z, dtype=np.float64)

    # Step 1: Apply the convolution (mean filter)
    kernel = np.ones((ker_size, ker_size)) / (ker_size ** 2)
    Z_m = convolve(Z, kernel, mode='reflect')  # Convolution to get smoothed surface

    # Step 2: Detrend by subtracting the smoothed surface (residual topography)
    Z_res = Z - Z_m

    # Step 3: Define the random sampling function for nlfilter-like behavior
    def sample_non_nan(values):
        # Sample random values, excluding NaN (represented by 0 in the matrix)
        valid_values = values[~np.isnan(values)]
        if len(valid_values) > 0:
            return np.median(random.sample(valid_values.tolist(), min(len(valid_values), n_samples)))
        else:
            return np.nan

    # Step 4: Apply the random sampling filter on the residual topography (Z_res)
    Z_window_random = generic_filter(Z_res, sample_non_nan, size=(big_window, big_window))

    return Z_window_random
")

# wrap it in an R function
detrend_surface <- function(Z, ker_size = 9L, big_window = 41L, n_samples = 1L) {
  
  m <- terra::as.matrix(Z, wide = TRUE) # Convert raster to matrix
  
  results <- py$detrend_surface(m, ker_size, big_window, n_samples) # Apply function
  
  r <- terra::rast(results, ext = ext(Z), crs = crs(Z)) # convert back to raster
  
  return(r)
}

##################
### SMALL TEST ###
##################

# library(terra)
# # Create raster
# r <- rast(nrows=100, ncols=100, xmin=0, xmax=100, ymin=0, ymax=100)
# values(r) <- runif(ncell(r), min=1, max=100)
# 
# # Introduce a missing data region
# r[40:60, 40:60] <- NA
# plot(r, main = "Raster with NA's")
# 
# # Make a mask layer for later
# A <- ifel(is.na(r), 0, r)
# # Fill in NA's with simple mean
# r <- ifel(is.na(r), mean(na.omit(values(r))), r)
# plot(r, main = "Raster with filled NA's")
# 
# # Apply the detrend to get a noise layer
# noise <- detrend_surface(r)
# 
# # Where there were originally NA's add noise to the filled values
# d <- ifel(A < 1, (r + noise), r)
# plot(d, main = "Raster with detrending")
