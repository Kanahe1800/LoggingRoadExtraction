# Theory — LoggingRoadExtraction

This document covers the theoretical background behind each stage of the pipeline. For setup and usage instructions, see [README.md](README.md).

---

## Stage 1 — Seed Generation (`createSeeds_v2.R`)

The seed generation stage produces a continuous road likelihood map (0–1) by combining three parallel terrain-derived feature streams. Each stream is independently normalized and fused into a single probability surface.

### 1. Slope Difference (`slope_dif`)

Raw slope captures local gradient, but roads do not simply have low slope — they have **locally low slope relative to their surroundings**. A large uniform low-pass filter (configurable via `SLOPE_SMOOTH_WINDOW`, default 51×51) estimates the regional slope trend. Subtracting this smoothed surface from the raw slope isolates pixels that are **anomalously flat** compared to their neighbourhood — a strong indicator of road cuts or fills.

$$slope_{dif} = slope - slope_{smooth}$$

Low (negative) values of `slope_dif` indicate road-like flatness. The window size is intentionally large to capture regional trend rather than local variation.

### 2. Multiscale Roughness (`roughness_ms`)

Roads are smooth. Surface roughness — estimated as the standard deviation of slope within a focal window — is low over road surfaces and high in forested or rocky terrain. To be robust across varying road widths and survey resolutions, roughness is computed at three spatial scales and averaged:

$$roughness_{ms} = \frac{\sigma_{small} + \sigma_{med} + \sigma_{large}}{3}$$

Window sizes are configurable via `ROUGHNESS_SD_SMALL`, `ROUGHNESS_SD_MED`, and `ROUGHNESS_SD_LARGE` (defaults: 5×5, 11×11, 21×21). Roads consistently exhibit **low roughness at all scales**, making this a multi-scale discriminator that is less sensitive to noise at any single scale.

### 3. Profile Curvature Edges (`edges`)

Profile curvature measures how slope changes along the direction of maximum gradient. Roads intersect natural terrain at sharp transitions (road shoulders and ditches), creating strong **linear curvature discontinuities**. These are detected in two steps:

1. A Sobel edge detector (configurable kernel size via `SOBEL_SIZE`, default 5×5) is applied to the profile curvature layer to compute gradient magnitude.
2. The resulting edge magnitudes are smoothed with a Gaussian kernel (`GAUSSIAN_SIZE` × `GAUSSIAN_SIZE`, σ = `GAUSSIAN_SIGMA`) to reduce noise and widen the signal around road boundaries.

$$edges = \mathcal{G}_{\sigma} \left( \| \nabla \, profcurv \| \right)$$

High edge magnitude indicates sharp terrain breaks consistent with road margins.

### Otsu Thresholding & Score Normalization

Each of the three feature layers is automatically thresholded using **Otsu's method**, which finds the optimal binary split by minimizing intra-class variance. Rather than using the binary result directly, the threshold serves as a **reference point** for soft normalization:

- For roughness and edges: pixels above the Otsu threshold are road-positive evidence.
- For slope difference: pixels *below* the Otsu threshold are road-positive evidence (anomalously flat).
- The margin from the threshold is scaled to [0, 1] by dividing by the maximum observed margin.

Otsu thresholding operates on an 8-bit rescaled histogram (0–255) for algorithm consistency, then back-transforms to real-value units.

This produces a continuous score per feature reflecting both **membership** (above/below threshold) and **confidence** (distance from threshold).

### Score Fusion

The three normalized scores are averaged into a single seed strength raster:

$$seeds = \frac{slope_{score} + rough_{score} + edges_{score}}{3}$$

Values near 1.0 indicate strong convergent evidence from all three terrain signals that a pixel lies on or near a road surface. The output is a continuous [0, 1] raster saved as `seeds_combined.tif`.

---

## Stage 2 — Raster to Vector (`RasterToVector.R`)

Converts the continuous seed surface into a sparse set of georeferenced candidate road nodes, snapped to true road-centre peaks, suitable for least-cost path routing.

### 1. Thresholding & Binary Mask

The seed surface is binarized at `PROB_THRESHOLD`. Pixels at or above this value are treated as candidate road pixels. An optional morphological **opening** (erosion via `focal min` followed by dilation via `focal max`, kernel size `OPENING_SIZE`) can remove isolated noise pixels before skeletonization.

### 2. Skeletonization (Zhang-Suen Thinning)

The binary road mask is reduced to a **1-pixel-wide skeleton** using the Zhang-Suen iterative thinning algorithm (`skeletonize.R`). The algorithm works by iteratively removing pixels from the boundary of the binary mask in two sub-iterations per pass, subject to three conditions:

1. The pixel has between 2 and 6 neighbours (prevents isolated pixel deletion and gap creation).
2. The sequence of alternating foreground/background neighbours around the pixel transitions exactly once (preserves connectivity).
3. Directional conditions differ between sub-iterations to prevent diagonal disconnection.

This collapses road-width blobs to their centrelines without breaking connectivity, ensuring downstream point sampling captures road position rather than road width.

### 3. Connected Component Filtering

The skeleton is labelled into 8-connected components via `terra::patches()`. Two filters are applied:

1. **Minimum blob size** (`MIN_SKEL_BLOB_PX`): removes isolated pixels caused by noise.
2. **Minimum component length** (`MIN_COMPONENT_LEN_M`): removes short fragments unlikely to represent actual roads. Length is estimated as `pixel count × pixel diagonal (m)`.

### 4. Global Focal Snapping

Raw skeleton pixels lie on the thinned centreline of the thresholded mask, which may not coincide with the highest-probability pixel (the true road centre) in the original `seeds` raster. Snapping corrects this in three steps:

1. A focal maximum filter is applied to the full `seeds` raster with a window of `(snap_radius × 2) + 1` pixels to identify local intensity peaks across the raster.
2. Local maxima are masked to a buffer around the skeleton (radius = `snap_radius + 1` pixels) to restrict the search space and keep RAM usage low.
3. Each skeleton point is matched to its nearest local maximum using `RANN::nn2()`, shifting it onto the highest-confidence road-centre pixel within the search window.

This produces `xy_snapped` — skeleton-derived points relocated to true road peaks.

### 5. Grid-Based Greedy Thinning

To produce a spatially uniform set of seed nodes, points are thinned using a **grid-cell approach**:

1. All snapped points are projected to UTM for metric coordinates.
2. Each point is assigned to a grid cell of size `MIN_SPACING_M × MIN_SPACING_M` metres using integer division of the coordinates.
3. Within each cell, only the point with the **highest seed score** is retained.

This is equivalent to greedy thinning with a regular spacing guarantee, but runs in O(n) time without a nearest-neighbour search — making it efficient on large skeletons with tens of thousands of points.

### CRS Handling

If the input raster is in geographic coordinates (lon/lat), it is automatically reprojected to the appropriate **UTM zone** (derived from the raster centroid) before any metric operations. The UTM zone is computed from the centroid longitude as `floor((lon + 180) / 6) + 1`, with hemisphere determined by latitude sign.

---

## Stage 2b — Interactive Point Editing (`PointEditingApp_v2.R`)

A Shiny application for manually reviewing and editing the candidate seed points before road extraction. It is not a processing stage in the algorithmic sense, but plays an important quality-control role.

### Design Rationale

Automated seed point generation cannot perfectly distinguish road pixels from other low-roughness, low-slope features (e.g. riverbeds, clearcuts, or flat rock outcrops). Manual review allows the operator to:

- **Delete** spurious points that fall outside known road corridors.
- **Add** points in road sections where the seed raster failed to detect the road (e.g. heavily vegetated or narrow segments).
- **Delete by AOI** — draw a polygon to bulk-delete all points within an area of interest, which is faster than clicking individual points in dense regions.

### Raster Overlay

The seed raster (`seeds_combined.tif`) is displayed as a background layer using `leafem::addRasterImage()`, reprojected to WGS84 for web rendering. This gives the operator visual context for where high road-likelihood areas are, allowing informed decisions about which points to keep or remove.

The raster is loaded **after** the Leaflet map initializes via `session$onFlushed`, which fires once the browser has received and rendered the initial map state. A `reactiveVal(FALSE)` flag (`map_ready`) gates the raster load to prevent `leafletProxy` being called before the map widget exists in the DOM — a common source of silent failures in Shiny-Leaflet applications.

### Map Interaction Model

- **Add mode**: map clicks create new `sf` point objects appended to the reactive point store (`rv`). Each point receives a unique incrementing integer ID.
- **Delete mode**: clicking a rendered circle marker fires `input$map_marker_click`, which filters the point out of `rv` by ID and removes the marker via `leafletProxy`.
- **AOI delete**: the Leaflet Draw toolbar captures polygon geometry as GeoJSON, which is converted to an `sfc` object and intersected with `rv` via `sf::st_intersects()`. All intersecting points are removed in a single operation.
- **Map view preservation**: the current zoom and centre are stored in a `map_view` reactive after every map interaction and restored via `setView()` after every proxy update, preventing the map from jumping back to the default view on each edit.

---

## Stage 3 — Road Extraction (`RoadExtraction.R`)

Connects candidate seed points via least-cost paths to reconstruct road centrelines as a vector layer.

### Cost Surface Construction

The seed raster is inverted to form a **resistance surface** — high seed values (likely roads) become low cost, and open terrain becomes high cost:

$$cost = 1 - seeds$$

NA cells are filled with the maximum observed cost value to ensure full raster coverage. A log transform is applied to compress the dynamic range and increase sensitivity in the low-cost (road) range:

$$cost_{log} = \left(\log(cost + \epsilon) - \min(\log(cost + \epsilon))\right) \times 100$$

The subtraction of the minimum and multiplication by 100 rescales the surface to start at 0, making threshold values interpretable as percentages of the dynamic range.

A **relative barrier** is applied using the maximum cost value — cells in the top 5% of the cost range are assigned a large penalty (`COST_BARRIER_VALUE`), and cells in the bottom 90% are further suppressed by a factor of 0.7 to increase contrast between road and non-road areas:

```r
cost_log <- ifel(cost_log > (cost_max * 0.95), COST_BARRIER_VALUE, cost_log)
cost_log <- ifel(cost_log < (cost_max * 0.90), cost_log * 0.7, cost_log)
```

### Conductance Transition Matrix

The resistance surface is inverted to **conductance** (higher = easier to traverse):

$$conductance = \frac{1}{cost_{log} + \epsilon}$$

The conductance surface is passed to `gdistance::transition()` with a **harmonic mean** transition function and 16-direction connectivity (`TR_DIRECTIONS`):

$$T_{ij} = \frac{2}{\frac{1}{c_i} + \frac{1}{c_j}}$$

The harmonic mean is preferred over the arithmetic mean because it is dominated by the *smaller* of the two cell conductances — meaning a single high-resistance cell on a path strongly penalizes that transition. This better enforces hard barriers and keeps paths confined to road corridors.

Geographic correction (`geoCorrection(type = "r")`) adjusts transition weights for the actual Euclidean distance between cell centres, accounting for the longer diagonal distance traversed in 8- and 16-direction grids.

### Least-Cost Path Chaining

The algorithm builds road chains greedily using a while loop over unused seed points:

1. Start a new chain from the first unused point.
2. Find all unused points within `MAX_DISTANCE_THRESHOLD` metres (computed as `MIN_SPACING_M × MAX_DISTANCE_THRESHOLD` from settings).
3. For each candidate in distance order, compute the least-cost path using `gdistance::shortestPath()` (Dijkstra's algorithm on the transition matrix).
4. Sample the resistance surface along each candidate path using `terra::extract()` and compute the mean cost.
5. Accept the candidate with the **lowest mean path cost** that is below `MEAN_COST_THRESHOLD`.
6. Mark the accepted point as used, append the path segment to the current chain, and repeat from the new endpoint.
7. When no viable connection is found within the search radius, close the current chain and start a new one from the next unused point.
8. Chains shorter than `MIN_SEGMENTS_PER_LINE` segments are discarded.

This greedy chaining approach approximates the branching structure of a logging road network without requiring a global graph solution. The lowest-cost selection — rather than nearest-neighbour — means the algorithm preferentially routes along road corridors even when a geometrically closer point would require crossing high-cost terrain.

### Viability Check

A connection is accepted only if the **mean resistance** sampled along the least-cost path is below `MEAN_COST_THRESHOLD`. This is evaluated on the original `cost_log` surface (resistance, not conductance), so lower values indicate the path stayed within road corridors. The check uses `terra::extract()` on the path geometry to sample all raster cells intersected by the line.