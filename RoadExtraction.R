library(terra)
library(sf)
library(gdistance)
library(raster)

dir.out <- "C:/Users/mimam/Dropbox/lider_data/output"

## Load data ##
seeds <- rast(file.path(dir.out, "seeds_combined.tif"))
points <- st_read(file.path(dir.out, "edited_points_2026-02-23.geojson"))

## Cost calculation
cost <- 1 - seeds
cost <- ifel(is.na(cost), max(values(cost, na.rm = TRUE)), cost)
cost_enhanced <- cost^2

# Add small value to avoid log(0)
cost_for_log <- cost + 0.001

# Log transform
cost_log <- log(cost_for_log)

## Check CRS
if (!is.na(st_crs(points)$input) && st_crs(points)$input != crs(seeds, proj = TRUE)) {
  points <- st_transform(points, crs = crs(seeds))
}

cost_raster <- raster(cost_enhanced)
tr <- transition(cost_raster, transitionFunction = mean, directions = 8)
tr <- geoCorrection(tr, type = "c", multpl = FALSE)

coords <- st_coordinates(points)

# Cost enhanced (squared)
plot(cost_enhanced, main = "Cost Enhanced (cost^2)\nUsed for Pathfinding")

# Log cost
plot(cost_log, main = "Cost Log (log(cost + 0.001))")


## Function to check if path is viable ##
is_path_viable <- function(lcp_sf, cost_terra, max_mean_cost = 0.5) {
  tryCatch({
    path_costs <- terra::extract(cost_terra, vect(lcp_sf), xy = TRUE)
    
    if (is.null(path_costs) || nrow(path_costs) == 0 || ncol(path_costs) < 2) {
      return(list(viable = FALSE, max_cost = Inf, mean_cost = Inf))
    }
    
    cost_values <- path_costs[[2]]
    cost_values <- cost_values[!is.na(cost_values)]
    
    if (length(cost_values) == 0) {
      return(list(viable = FALSE, max_cost = Inf, mean_cost = Inf))
    }
    
    max_cost <- max(cost_values)
    mean_cost <- mean(cost_values)
    
    viable <- !is.na(mean_cost) && !is.infinite(mean_cost) && mean_cost <= max_mean_cost
    
    return(list(viable = viable, max_cost = max_cost, mean_cost = mean_cost))
  }, error = function(e) {
    return(list(viable = FALSE, max_cost = Inf, mean_cost = Inf))
  })
}

## PARAMETERS ##
mean_cost_threshold <- 0.85
max_distance_threshold <- 300  # Maximum Euclidean distance to consider
min_segments_per_line <- 3     # Minimum segments to keep a line

print(paste("Total points:", nrow(coords)))

## BUILD CONNECTED ROAD CHAINS ##
print("\n========== BUILDING CONNECTED ROAD CHAINS ==========")

# Track which points have been used
used_points <- rep(FALSE, nrow(coords))

# Store all road lines
all_road_lines <- list()
line_counter <- 0

while (sum(!used_points) > 0) {
  # Start a new road line from an unused point
  available_points <- which(!used_points)
  if (length(available_points) == 0) break
  
  # Start from a random unused point (or first available)
  current_point_idx <- available_points[1]
  used_points[current_point_idx] <- TRUE
  
  current_line_segments <- list()
  segment_counter <- 0
  
  print(paste("\n--- Starting new road line", line_counter + 1, "from point", current_point_idx, "---"))
  
  # Build chain by repeatedly finding nearest viable neighbor
  while (TRUE) {
    # Find nearest unused point
    current_coords <- coords[current_point_idx, , drop = FALSE]
    
    # Calculate distances to all unused points
    unused_indices <- which(!used_points)
    if (length(unused_indices) == 0) break
    
    distances <- sqrt((coords[unused_indices, 1] - current_coords[1])^2 + 
                        (coords[unused_indices, 2] - current_coords[2])^2)
    
    # Sort by distance
    sorted_indices <- order(distances)
    
    # Try to connect to nearest points (try top 5)
    connected <- FALSE
    for (try_idx in sorted_indices[1:min(5, length(sorted_indices))]) {
      next_point_idx <- unused_indices[try_idx]
      next_dist <- distances[try_idx]
      
      # Skip if too far
      if (next_dist > max_distance_threshold) next
      
      # Try to create path
      tryCatch({
        lcp <- shortestPath(tr, 
                            origin = coords[current_point_idx, , drop = FALSE], 
                            goal = coords[next_point_idx, , drop = FALSE], 
                            output = "SpatialLines")
        
        if (!is.null(lcp)) {
          lcp_sf <- st_as_sf(lcp)
          st_crs(lcp_sf) <- st_crs(points)
          
          path_check <- is_path_viable(lcp_sf, cost_enhanced, max_mean_cost = mean_cost_threshold)
          
          if (path_check$viable) {
            # Accept this connection
            segment_counter <- segment_counter + 1
            current_line_segments[[segment_counter]] <- lcp
            used_points[next_point_idx] <- TRUE
            
            print(paste("  Segment", segment_counter, ": Point", current_point_idx, "->", 
                        next_point_idx, "| Dist:", round(next_dist, 1), 
                        "| Mean cost:", round(path_check$mean_cost, 3)))
            
            # Move to next point
            current_point_idx <- next_point_idx
            connected <- TRUE
            break
          }
        }
      }, error = function(e) {
        # Try next candidate
      })
    }
    
    # If no viable connection found, end this line
    if (!connected) {
      print(paste("  No more viable connections - ending line"))
      break
    }
  }
  
  # Save line if it meets minimum segment requirement
  if (segment_counter >= min_segments_per_line) {
    line_counter <- line_counter + 1
    all_road_lines[[line_counter]] <- current_line_segments
    print(paste("*** Saved road line", line_counter, "with", segment_counter, "segments ***"))
  } else {
    print(paste("*** Discarded line with only", segment_counter, "segments (< ", min_segments_per_line, ") ***"))
  }
}

print(paste("\n=== SUMMARY ==="))
print(paste("Total road lines created:", line_counter))
print(paste("Used points:", sum(used_points), "/", nrow(coords)))

## EXPORT RESULTS ##
if (line_counter > 0) {
  # Flatten all segments
  all_path_segments <- list()
  segment_idx <- 0
  
  for (line_id in 1:length(all_road_lines)) {
    for (seg in all_road_lines[[line_id]]) {
      segment_idx <- segment_idx + 1
      all_path_segments[[segment_idx]] <- seg
    }
  }
  
  all_paths_sp <- do.call(rbind, all_path_segments)
  all_paths_sf <- st_as_sf(all_paths_sp)
  st_crs(all_paths_sf) <- st_crs(points)
  
  # Add line ID attribute
  all_paths_sf$line_id <- NA
  segment_idx <- 1
  for (line_id in 1:length(all_road_lines)) {
    n_segs <- length(all_road_lines[[line_id]])
    all_paths_sf$line_id[segment_idx:(segment_idx + n_segs - 1)] <- line_id
    segment_idx <- segment_idx + n_segs
  }
  
  st_write(all_paths_sf, 
           file.path(dir.out, "logging_road_lcp_connected_chains.geojson"),
           delete_dsn = TRUE)
  
  print("\nConnected road chains exported successfully!")
  
  # Plot results with different colors per line
  colors <- rainbow(line_counter)
  
  plot(cost_enhanced, main = paste("Connected Road Chains\n", 
                                   line_counter, "lines,", 
                                   segment_idx - 1, "segments"))
  
  for (line_id in 1:line_counter) {
    line_segs <- all_paths_sf[all_paths_sf$line_id == line_id, ]
    plot(st_geometry(line_segs), add = TRUE, col = "red", lwd = 3)
  }
  
  # Plot used vs unused points
  used_pts <- points[used_points, ]
  unused_pts <- points[!used_points, ]
  
  # plot(st_geometry(used_pts), add = TRUE, col = "blue", pch = 16, cex = 1)
  # if (nrow(unused_pts) > 0) {
  #   plot(st_geometry(unused_pts), add = TRUE, col = "orange", pch = 17, cex = 1.5)
  # }
  
  # Print line statistics
  print("\n=== LINE STATISTICS ===")
  for (line_id in 1:line_counter) {
    n_segs <- sum(all_paths_sf$line_id == line_id)
    print(paste("Line", line_id, ":", n_segs, "segments"))
  }
  
} else {
  print("No road lines were created!")
}