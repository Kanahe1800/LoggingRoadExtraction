# Limitations — LoggingRoadExtraction

This document describes the known limitations of the current pipeline. Some are fundamental to the underlying methods, while others are active areas of development. For theoretical background, see [THEORY.md](THEORY.md).

---

## Stage 1 — Seed Generation (`createSeeds_v2.R`)

### Flat terrain performance
The pipeline relies heavily on **slope difference** (`slope_dif`) to identify anomalously flat areas relative to the surrounding landscape. In naturally flat terrain — such as valley bottoms, plateaus, agricultural areas, or lakeshores — the regional slope trend is already low, so roads produce little detectable contrast against the background. This causes the slope difference stream to contribute weak or noisy signal, reducing overall seed quality in flat areas.

### Watershed and stream channel confusion
Small **stream channels and riverbeds** share many of the same terrain signatures as logging roads: they are locally flat relative to their surroundings, exhibit low multiscale roughness, and create sharp curvature transitions at their banks. As a result, the seed raster frequently assigns high road-likelihood scores to stream corridors. 

> **Status:** Masking of water features (using flow accumulation, wetness index, or external water body layers) is currently in development and not yet integrated into the pipeline.

### Clearcut and harvested area confusion
Recently harvested areas and clearcuts are **smooth, flat, and have sharp terrain edges** at their boundaries — all of which score highly across all three feature streams. These areas can produce broad high-seed zones that are difficult to distinguish from road networks without additional contextual layers (e.g. cutblock polygons or NDVI from optical imagery).

### Sensitivity to DEM resolution and quality
The pipeline was developed and tested on LiDAR-derived DEMs. Performance on **coarser or photogrammetric DEMs** (e.g. from drone SfM or SRTM) may be degraded because:
- Road widths may be sub-pixel at coarse resolutions, making skeletonization unreliable.
- Photogrammetric DEMs can have surface noise that mimics road roughness signatures.
- Canopy penetration is limited in optical DEMs, which may suppress road signals under forest cover.

### Fixed kernel sizes
The slope smoothing window (`SLOPE_SMOOTH_WINDOW`), roughness scales, and Gaussian kernel are set as fixed pixel sizes in `settings.env`. These were tuned for a specific DEM resolution. If the **pixel size of your DEM differs significantly** from the development dataset, these parameters should be scaled proportionally — the pipeline does not currently do this automatically.

### Profile curvature sensitivity to noise
The Sobel edge detector applied to profile curvature is sensitive to **high-frequency noise** in the DEM. In DEMs with significant acquisition noise or artefacts (e.g. flightline striping in LiDAR), this can generate spurious edge features that inflate seed scores in non-road areas.

### Equal weighting of feature streams
The three feature streams (slope, roughness, edges) are combined with **equal weights** (1/3 each). In some terrain types, one stream may be far more informative than the others — for example, edges may dominate in steep terrain while slope difference is most useful in gentle terrain. There is currently no adaptive or learned weighting scheme.

---

## Stage 2 — Raster to Vector (`RasterToVector.R`)

### Skeletonization breaks at road junctions
The Zhang-Suen thinning algorithm is designed to produce **1-pixel-wide connected skeletons**, but it can fail at road junctions and T-intersections, sometimes breaking the skeleton into disconnected fragments or creating spurious branches. This leads to missing seed points at junction areas.

### Snapping may pull points off-road
The focal snapping step relocates skeleton points to the nearest local maximum within `snap_radius` pixels. If the **local maximum is not on the road** (e.g. it falls on a nearby ridge or stream bank that scored highly), points can be incorrectly displaced away from the true road centreline.

### Grid thinning is not spatially adaptive
The grid-based thinning assigns one point per `MIN_SPACING_M × MIN_SPACING_M` grid cell, keeping the highest-scoring point. This means **winding or curved roads** that cross multiple grid cells may be well-represented, but **straight roads** along a single grid column or row may be undersampled, producing fewer seed points than needed for the LCP chaining stage.

### Hard-coded UTM zone
The current implementation defaults to **EPSG:32610** (UTM Zone 10N) for metric operations during thinning. This is correct for coastal British Columbia but will produce incorrect distances for DEMs outside this zone. The `to_utm_epsg()` function computes the correct zone dynamically, but the thinning step still uses the hardcoded value.

> **Status:** This should be updated to use the dynamically computed UTM zone throughout.

### Memory usage on large DEMs
The focal maximum filter applied during snapping operates on the **full raster extent**. For very large DEMs, this can require substantial RAM. The current workaround (masking maxima to the skeleton buffer) reduces but does not eliminate this issue.

---

## Stage 2b — Interactive Point Editing (`PointEditingApp_v2.R`)

### Manual editing is operator-dependent
The quality of the final road extraction depends heavily on **operator judgement** during the editing step. There are no automated checks for whether added points are placed on roads or deleted points were actually spurious. Results are not reproducible across different operators without editing guidelines.

### No undo functionality
Point edits cannot be undone individually. The only recovery option is **Reset to uploaded**, which discards all edits made since the last file load. There is no edit history or incremental undo stack.

### Performance degrades with large point sets
The app renders all points as individual `CircleMarker` objects in the browser. With **thousands of points**, rendering and interaction can become slow. For large study areas, it is recommended to increase `MIN_SPACING_M` to reduce point density before launching the app.

### Raster overlay is slow for large rasters
`leafem::addRasterImage()` transfers the full raster to the browser as an image. For **large or high-resolution rasters**, this can cause a long load time or browser memory issues. Downsampling with `aggregate(cost_raster, fact = 5)` (commented out in the script) can help at the cost of visual detail.

### No coordinate display
The app does not show the coordinates of the cursor or clicked points. This makes it difficult to precisely place points at known road locations without cross-referencing with an external GIS application.

---

## Stage 3 — Road Extraction (`RoadExtraction.R`)

### Greedy chaining does not guarantee global optimality
The road chain building algorithm is **greedy** — it always selects the locally lowest-cost connection from the current endpoint. This means it can make early connections that prevent better global routing later. A true global graph solution (e.g. minimum spanning tree or network flow) would produce better results but is computationally more expensive.

### Sensitivity to seed point placement
The LCP algorithm connects whatever points are provided. If seed points are **poorly placed** (too sparse, incorrectly located, or missing at key junctions), the extracted road lines will be incomplete or incorrectly routed regardless of the cost surface quality.

### No junction detection
The chaining algorithm builds **linear chains** — it does not detect or reconstruct road junctions, T-intersections, or branching points. Each chain starts from an unused point and grows linearly until no viable connection is found. As a result, complex road networks with many intersections are often split into disconnected segments.

### Road direction is not considered
Least-cost paths are computed **symmetrically** — the cost from point A to point B is the same as from B to A. In reality, road construction follows terrain in ways that may be directionally asymmetric (e.g. preferring to follow contours in one direction). An anisotropic cost surface would better capture this, but is not currently implemented.

### Performance scales poorly with point count
For each candidate connection, the algorithm calls `shortestPath()`, which runs Dijkstra's algorithm on the full transition matrix. With many seed points and a large `MAX_DISTANCE_THRESHOLD`, the number of LCP computations grows quadratically. **Large study areas with dense seed points can take hours to process.**

### Cost surface is relative, not absolute
The cost thresholds (`MEAN_COST_THRESHOLD`, `COST_BARRIER_VALUE`) are defined as absolute values on a normalized scale, but the normalization is relative to each individual raster's dynamic range. This means **thresholds are not transferable between datasets** — values tuned for one study area may need significant adjustment for another.

### 16-direction transition can still produce diagonal artefacts
While 16-direction connectivity (`TR_DIRECTIONS = 16`) substantially reduces the staircase artefacts seen with 8-direction grids, **some diagonal stepping** may still be visible in paths crossing open terrain, particularly where the cost surface is uniform and the algorithm has no directional preference.

---

## General Limitations

### Memory intensity — not suitable for large DEMs
The pipeline is **memory-intensive at every stage** and is not currently designed to handle large DEM extents without manual tiling or downsampling:

- **Stage 1** loads the full DEM and computes multiple focal filter operations (slope, roughness at three scales, profile curvature, Sobel edges, Gaussian smoothing) entirely in RAM. A 51×51 focal window on a large raster can require tens of gigabytes of memory. Intermediate objects are cleared with `rm()` + `gc()` between steps, but peak RAM usage can still exceed available memory on typical workstations.
- **Stage 2** holds the full `seeds` raster, skeleton matrix, focal maximum raster, and point datasets simultaneously in memory during the snapping step. For large rasters this can cause R to crash or swap heavily to disk.
- **Stage 3** builds a `gdistance` transition matrix over the full raster extent. The transition matrix is a sparse matrix whose size scales with the number of raster cells — for large DEMs this object alone can exhaust available RAM before any paths are computed.

**Recommended workarounds (not yet automated):**
- Clip the input DEM to the study area of interest before running the pipeline.
- Tile the DEM into overlapping blocks, process each tile independently, and merge outputs.
- Downsample the DEM to a coarser resolution if sub-metre precision is not required.
- Increase `MIN_SPACING_M` to reduce the number of seed points passed to Stage 3.

> **Status:** Tiled processing and out-of-memory raster handling are not currently implemented. The pipeline is best suited to DEM extents of a few square kilometres at LiDAR resolution (≤1m pixel).

### No validation against ground truth
The pipeline has no built-in accuracy assessment. There is currently no automated comparison against reference road datasets (e.g. forest tenure road layers or GPS-surveyed centrelines). Accuracy must be assessed manually by the operator using visual inspection or external GIS tools.

### Single-date DEM only
The pipeline uses a single LiDAR acquisition. It cannot detect **road abandonment, revegetation, or new road construction** that occurred after the survey date. Time-series analysis would require multiple DEMs or integration with optical satellite imagery.

### No canopy height model integration
Roads under **dense forest canopy** may be partially obscured in the DEM, reducing slope and roughness contrast. Integrating a Canopy Height Model (CHM) to mask or down-weight high-canopy areas could improve detection in closed-canopy environments, but this is not currently implemented.

### Projection assumed to be metric
All distance calculations assume the raster is projected in a metric CRS (e.g. UTM). If the input DEM is in geographic coordinates (degrees), pixel-based distance estimates will be incorrect. The pipeline includes an automatic reprojection step in `RasterToVector.R`, but `RoadExtraction.R` does not currently verify CRS before computing Euclidean distances between seed points.