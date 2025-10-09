## INPAINT_NANS FUNCTION WRITTEN IN PYTHON ##
## translated method 0/1 created by: John D'Errico (2025). inpaint_nans (https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans), MATLAB Central File Exchange. Retrieved March 11, 2025.

################
### COMMENTS ###
################

#* in the future I may need to set this up so that the function is in a .py
#* file which acts as a source for an r function.
#*
#* need to look into further optimizing the python function to reduce mem. usage
#* 
#* using float64 also seems overkill and I believe introduces zeros after the
#* decimal which might needlessly increase the size of the output layer
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

# install the dependant python package(s)
# py_install("scipy") # if required

###########################
### DEFINE THE FUNCTION ###
###########################

# define in python
py_run_string("
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

def inpaint_nans(A, method=0):
    if method not in [0, 1]:
        raise ValueError('Only methods 0 and 1 are supported.')

    A = np.array(A, dtype=np.float64)
    nan_mask = np.isnan(A)
    if not np.any(nan_mask):
        return A

    rows, cols = A.shape
    idx_matrix = np.arange(rows * cols).reshape(rows, cols)
    A_flat = A.flatten()

    I, J, V = [], [], []  # Sparse matrix components
    b = np.zeros(rows * cols)

    for r in range(rows):
        for c in range(cols):
            idx = idx_matrix[r, c]
            if not nan_mask[r, c]:
                I.append(idx)
                J.append(idx)
                V.append(1)
                b[idx] = A_flat[idx]
            else:
                neighbors = []
                for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        neighbors.append(idx_matrix[nr, nc])

                if method == 0:
                    weight = -1 / len(neighbors)
                    for n_idx in neighbors:
                        I.append(idx)
                        J.append(n_idx)
                        V.append(weight)
                    I.append(idx)
                    J.append(idx)
                    V.append(1)

                elif method == 1:
                    weight = -1 / (4 + len(neighbors))
                    for n_idx in neighbors:
                        I.append(idx)
                        J.append(n_idx)
                        V.append(weight)
                    I.append(idx)
                    J.append(idx)
                    V.append(1)

    A_sparse = lil_matrix((rows * cols, rows * cols))
    A_sparse[I, J] = V
    x = spsolve(A_sparse.tocsr(), b)

    return x.reshape(rows, cols)
")

# wrap it in an R function
inpaint_nans <- function(x, method = 0) {
  
  m <- terra::as.matrix(x, wide = TRUE) # Convert to a matrix
  
  # Apply the python function
  results <- py$inpaint_nans(m, method)
  
  r <- terra::rast(results, ext = ext(x), crs = crs(x)) # Convert back to raster
  
  return(r)
}

##################
### SMALL TEST ###
##################

# library(terra)
# # Create raster
# r <- rast(nrows=100, ncols=100, xmin=0, xmax=100, ymin=0, ymax=100)
# values(r) <- runif(ncell(r), min=0, max=100)
# 
# # Introduce a missing data region
# r[40:60, 40:60] <- NA
# plot(r, main = "Raster with NA's")
# 
# r_filled <- inpaint_nans(r)
# plot(r_filled, main = "Raster with filled NA's")
