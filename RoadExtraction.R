library(terra)
library(sf)
library(gdistance)
library(raster)

source("EnvLoader.R")

## Load data ##
seeds <- rast(seeds_file)
points <- st_read(edited_points_file)

## ── Cost surface ──────────────────────────────────────────────────────────
cost          <- 1 - seeds
cost          <- ifel(is.na(cost), max(values(cost, na.rm = TRUE)), cost)
cost_log <- log(cost + 0.001)
cost_log <- (cost_log - min(values(cost_log, na.rm = TRUE))) * 100
cost_max = max(values(cost_log, na.rm = TRUE))
print(cost_max)
cost_log <- ifel(cost_log > (cost_max - cost_max * 0.05), cost_barrier_value, cost_log)
cost_log <- ifel(cost_log < (cost_max - cost_max * 0.1), cost_log * 0.7, cost_log)

conductance <- 1 / (cost_log + 0.001)   # invert: low cost → high conductance
cost_raster <- raster(conductance)
plot(cost_log, main = "Cost Log Inversed (log(cost + 0.001))")

## ── Transition matrix ─────────────────────────────────────────────────────
tr <- transition(cost_raster, 
                 transitionFunction = function(x) 2 / (1/x[1] + 1/x[2]),
                 directions = tr_directions)
tr <- geoCorrection(tr, type = "r", multpl = FALSE)

## ── Align CRS ─────────────────────────────────────────────────────────────
if (!is.na(st_crs(points)$input) &&
    st_crs(points)$input != crs(seeds, proj = TRUE)) {
  points <- st_transform(points, crs = crs(seeds))
}

coords <- st_coordinates(points)
n      <- nrow(coords)
print(paste("Total points:", n))

## ── PARAMETERS ────────────────────────────────────────────────────────────
max_distance_threshold <- min_spacing_m * max_distance_threshold    # Search radius (m) — try all unused points within this
mean_cost_threshold    <- mean_cost_threshold   # Max mean cost^2 along path to accept connection
min_segments_per_line  <- min_segments_per_line      # Discard chains shorter than this

## ── Helper: viability check ───────────────────────────────────────────────
is_path_viable <- function(lcp_sf, cost_terra, max_mean_cost) {
  tryCatch({
    path_costs  <- terra::extract(cost_terra, vect(lcp_sf), xy = TRUE)
    cost_values <- path_costs[[2]]
    cost_values <- cost_values[!is.na(cost_values)]
    if (length(cost_values) == 0)
      return(list(viable = FALSE, mean_cost = Inf))
    mean_cost <- mean(cost_values)
    print(mean_cost)
    list(viable    = !is.na(mean_cost) && mean_cost <= max_mean_cost,
         mean_cost = mean_cost)
  }, error = function(e) list(viable = FALSE, mean_cost = Inf))
}

## ── BUILD CONNECTED ROAD CHAINS ───────────────────────────────────────────
print("\n========== BUILDING CONNECTED ROAD CHAINS ==========")

used_points    <- rep(FALSE, n)
all_road_lines <- list()
line_counter   <- 0

while (sum(!used_points) > 0) {
  
  available_points    <- which(!used_points)
  current_point_idx   <- available_points[1]
  used_points[current_point_idx] <- TRUE
  
  current_line_segments <- list()
  segment_counter       <- 0
  
  print(paste("\n--- Starting new road line", line_counter + 1,
              "from point", current_point_idx, "---"))
  
  repeat {
    current_coords <- coords[current_point_idx, , drop = FALSE]
    unused_indices <- which(!used_points)
    if (length(unused_indices) == 0) break
    
    ## Distances to all unused points
    distances <- sqrt((coords[unused_indices, 1] - current_coords[1])^2 +
                        (coords[unused_indices, 2] - current_coords[2])^2)
    
    ## ── KEY CHANGE: consider ALL unused points within the distance limit ──
    ## Sort by distance (nearest first), but don't cap at top-5
    candidates <- unused_indices[order(distances)]
    cand_dists <- distances[order(distances)]
    in_range   <- cand_dists <= max_distance_threshold
    
    candidates <- candidates[in_range]
    cand_dists <- cand_dists[in_range]
    
    if (length(candidates) == 0) {
      print("  No unused points within distance limit — ending line")
      break
    }
    
    ## Try every candidate in distance order; pick the first viable one
    ## (viable = LCP mean cost below threshold)
    best_idx      <- NA
    best_lcp      <- NULL
    best_mean_cost <- Inf
    best_dist     <- Inf
    
    for (k in seq_along(candidates)) {
      next_point_idx <- candidates[k]
      next_dist      <- cand_dists[k]
      
      tryCatch({
        lcp <- shortestPath(tr,
                            origin = current_coords,
                            goal   = coords[next_point_idx, , drop = FALSE],
                            output = "SpatialLines")
        if (is.null(lcp)) return()
        
        lcp_sf <- st_as_sf(lcp)
        st_crs(lcp_sf) <- st_crs(points)
        
        chk <- is_path_viable(lcp_sf, cost_log,
                              max_mean_cost = mean_cost_threshold)
        
        ## LOWEST-COST: inside the loop, no break
        if (chk$viable && chk$mean_cost < best_mean_cost) {
          best_idx       <- next_point_idx
          best_lcp       <- lcp_sf
          best_mean_cost <- chk$mean_cost
          best_dist      <- next_dist
        }
      }, error = function(e) {})
    }
    
    # ── If you prefer the LOWEST-COST viable path instead of nearest,
    #    remove the `break` above and replace with:

      if (chk$viable && chk$mean_cost < best_mean_cost) {
        best_idx       <- next_point_idx
        best_lcp       <- lcp_sf
        best_mean_cost <- chk$mean_cost
        best_dist      <- next_dist
      }

    ## (then accept best_idx after the loop)
    
    if (is.na(best_idx)) {
      print("  No viable connection among all candidates — ending line")
      break
    }
    
    ## Accept connection
    segment_counter <- segment_counter + 1
    current_line_segments[[segment_counter]] <- best_lcp
    used_points[best_idx] <- TRUE
    
    print(paste("  Segment", segment_counter,
                ": Point", current_point_idx, "->", best_idx,
                "| Dist:", round(best_dist, 1),
                "| Mean cost:", round(best_mean_cost, 3)))
    
    current_point_idx <- best_idx
  }
  
  ## Save chain if long enough
  if (segment_counter >= min_segments_per_line) {
    line_counter <- line_counter + 1
    all_road_lines[[line_counter]] <- current_line_segments
    print(paste("*** Saved road line", line_counter,
                "with", segment_counter, "segments ***"))
  } else {
    print(paste("*** Discarded line with", segment_counter,
                "segments (< min", min_segments_per_line, ") ***"))
  }
}

print(paste("\n=== SUMMARY ==="))
print(paste("Total road lines created:", line_counter))
print(paste("Points used:", sum(used_points), "/", n))
print(paste("Points unused (isolated):", sum(!used_points)))

## ── Export ────────────────────────────────────────────────────────────────
if (line_counter > 0) {
  
  all_path_segments <- list()
  segment_idx       <- 0

  for (line_id in seq_len(line_counter)) {
    for (seg in all_road_lines[[line_id]]) {
      segment_idx <- segment_idx + 1
      all_path_segments[[segment_idx]] <- seg
    }
  }

  all_paths_sf <- do.call(rbind, all_path_segments)
  st_crs(all_paths_sf) <- st_crs(points)

  ## Add line_id attribute
  all_paths_sf$line_id <- NA_integer_
  seg_start <- 1
  for (line_id in seq_len(line_counter)) {
    n_segs <- length(all_road_lines[[line_id]])
    all_paths_sf$line_id[seg_start:(seg_start + n_segs - 1)] <- line_id
    seg_start <- seg_start + n_segs
  }

  if (file.exists(final_road_path_file)) file.remove(final_road_path_file)
  st_write(all_paths_sf, final_road_path_file, driver = "GeoJSON", delete_dsn = TRUE)
  print(paste("\nExported:", final_road_path_file))
  
  ## ── Plot ────────────────────────────────────────────────────────────────
  plot(cost_log,
       main = paste("Connected Road Chains v2\n",
                    line_counter, "lines |",
                    segment_idx, "segments |",
                    sum(!used_points), "isolated points"))
  
  for (line_id in seq_len(line_counter)) {
    line_segs <- all_paths_sf[all_paths_sf$line_id == line_id, ]
    plot(st_geometry(line_segs), add = TRUE,
         col = "red", lwd = 2)
  }
  
  plot(st_geometry(points[used_points,  ]), add = TRUE,
       col = "orange",   pch = 16, cex = 0.7)
  # plot(st_geometry(points[!used_points, ]), add = TRUE,
  #      col = "orange", pch = 17, cex = 1.2)
  
  ## Line statistics
  # print("\n=== LINE STATISTICS ===")
  # for (line_id in seq_len(line_counter)) {
  #   for (seg in all_road_lines[[line_id]]) {
  #     seg$line_id    <- line_id        # tag here
  #     segment_idx    <- segment_idx + 1
  #     all_path_segments[[segment_idx]] <- seg
  #   }
  # }
  
} else {
  print("No road lines were created! Try raising max_distance_threshold or mean_cost_threshold.")
}