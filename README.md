# LoggingRoadExtraction


## Stage 1 - Seed Generation (`createSeeds.R`)

A LiDAR-based road detection pipeline that generates a continuous **seed likelihood map** (0–1) from a Digital Elevation Model (DEM). The seed layer identifies candidate road pixels by combining terrain-derived evidence across three parallel feature streams, each independently normalized and fused into a single probability surface.

### Workflow
![createSeeds workflow diagram](images/createSeeds.png)


###  Theory

#### 1. Slope Difference (`slope_dif`)
Raw slope captures local gradient, but roads do not simply have low slope — they have **locally low slope relative to their surroundings**. A large uniform low-pass filter (51×51 kernel) estimates the regional slope trend. Subtracting this smoothed surface from the raw slope isolates pixels that are **anomalously flat** compared to their neighbourhood, which is a strong indicator of road cuts or fills.

$$\text{slope\_dif} = \text{slope} - \text{slope\_smooth}$$

Low (negative) values of `slope_dif` indicate road-like flatness.

#### 2. Multiscale Roughness (`roughness_ms`)

Roads are smooth. Surface roughness — estimated as the standard deviation of slope within a focal window — is low over road surfaces and high in forested or rocky terrain. To be robust across varying road widths and contexts, roughness is computed at three spatial scales (5×5, 11×11, 21×21 pixels) and averaged:

$$\text{roughness\_ms} = \frac{\sigma_5 + \sigma_{11} + \sigma_{21}}{3}$$

Roads consistently exhibit **low roughness at all scales**, making this a multi-scale discriminator that is less sensitive to noise at any single scale.

#### 3. Profile Curvature Edges (`edges`)

Profile curvature measures how the slope is changing along the direction of maximum gradient. Roads intersect natural terrain at sharp transitions (road shoulders/ditches), creating strong **linear curvature discontinuities**. These are detected by applying a 5×5 Sobel edge detector to the profile curvature layer. The resulting edge magnitudes are smoothed with a 21×21 Gaussian kernel (σ=3) to reduce noise and widen the signal around road boundaries:

$$\text{edges} = \mathcal{G}_{\sigma=3} \left( \| \nabla \text{profcurv} \| \right)$$

High edge magnitude indicates sharp terrain breaks consistent with road margins.

#### Otsu Thresholding & Score Normalization

Each of the three feature layers is automatically thresholded using **Otsu's method**, which finds the optimal binary split by minimizing intra-class variance. Rather than using the binary result directly, the threshold is used as a **reference point** for soft normalization:

- Values above (or below, for slope) the Otsu threshold are treated as road-positive evidence.
- The margin from the threshold is scaled to [0, 1] by dividing by the maximum observed margin.

This produces a continuous *score* per feature that reflects both **membership** (above/below threshold) and **confidence** (how far above/below).

#### Score Fusion

The three normalized scores are averaged into a single seed strength raster:

$$\text{seeds} = \frac{\text{slope\_score} + \text{rough\_score} + \text{edges\_score}}{3}$$

Values near 1.0 indicate strong, convergent evidence from all three terrain signals that a pixel lies on or near a road surface. This continuous output is suitable for use as a cost surface or seed layer in subsequent least-cost path analysis.

### Dependencies
| Package | Purpose |
|---|---|
| `terra` | Raster I/O, focal filtering, terrain analysis |
| `whitebox` | Profile curvature via WhiteboxTools |
| `autothresholdr` | Otsu automatic thresholding |

---
### Directory Structure

```
project/
├── input/
│   └── your_LiDAR_data.tif          # Input LiDAR DEM
├── output/
│   ├── profcurv_cc.tif        # Intermediate: profile curvature
│   └── seeds_combined.tif     # Final output: seed likelihood map
└── tmp/                       # Temporary working files
```

---

### Configuration for createSeeds.R

Edit the top of `createSeeds.R` to set your paths and input file:

```r
dir.in   <- "path/to/input"
dir.out  <- "path/to/output"
dir.tmp  <- "path/to/tmp"
file_name <- "your_dem.tif"
```

---
### Output from createSeeds.R

`seeds_combined.tif` — a single-band raster co-registered with the input DEM, with pixel values in [0, 1] representing road likelihood. Higher values indicate stronger convergent evidence of road presence across slope, roughness, and curvature edge features.


## Stage 2 — Raster to Vector (`RasterToVector.R`)

Converts the `seeds` probability raster into a **thinned set of georeferenced seed points** (GeoJSON) suitable for use as inputs to least-cost path analysis or manual digitizing. The core function is `seeds_to_points()`.

### Workflow
![RasterToVector workflow diagram](images/rasterToVector.png)

### Theory

#### Thresholding & Binary Mask

The continuous seed surface is binarized at a user-defined `prob_threshold` (default 0.3). Pixels at or above this value are treated as candidate road pixels. An optional morphological **opening** (erosion followed by dilation) can be applied to remove isolated noise pixels before skeletonization.

#### Skeletonization (Zhang-Suen Thinning)

The binary road mask is reduced to a **1-pixel-wide skeleton** using the Zhang-Suen iterative thinning algorithm (sourced from `skeletonize.R`). This collapses road polygons to their centrelines, ensuring that downstream point sampling captures road position rather than road width.

#### Connected Component Filtering

The skeleton is labelled into 8-connected components. Two filters are applied sequentially:

1. **Minimum blob size** (`min_skel_blob_px`): removes isolated skeleton pixels likely caused by noise.
2. **Minimum component length** (`min_component_len_m`): removes short skeleton fragments unlikely to represent actual roads. Component length is estimated as pixel count × pixel diagonal (in metres).

#### Score-Weighted Point Thinning

Skeleton pixels are converted to point geometries and each point is assigned the corresponding seed score from the original raster. Points are then **greedily thinned** to a minimum spacing of `min_spacing_m` metres: points are processed in descending score order, and any neighbour within the spacing radius is suppressed. This produces a spatially distributed set of high-confidence seed points using a fixed-radius nearest-neighbour index (`dbscan::frNN`) for efficiency.

#### CRS Handling

If the input raster is in geographic coordinates (lon/lat), it is automatically reprojected to the appropriate **UTM zone** (derived from the raster centroid) before any metric operations. Output points are re-projected back to WGS84 (EPSG:4326) before saving.

### Dependencies

| Package | Purpose |
|---|---|
| `terra` | Raster operations, coordinate extraction |
| `sf` | Vector I/O, CRS transforms, GeoJSON output |
| `dbscan` | Fixed-radius nearest-neighbour search for point thinning |

`skeletonize.R` (local source) must be present and export `thin_zhang_suen()`.

### Configuration

Key parameters in the `seeds_to_points()` call:

| Parameter | Default | Description |
|---|---|---|
| `prob_threshold` | `0.25` | Minimum seed score to include a pixel |
| `do_opening` | `FALSE` | Apply morphological opening before skeletonization |
| `opening_size` | `3` | Kernel size for opening (pixels) |
| `min_skel_blob_px` | `1` | Minimum skeleton fragment size (pixels) |
| `min_component_len_m` | `1.0` | Minimum skeleton component length (metres) |
| `min_spacing_m` | `8.0` | Minimum spacing between output points (metres) |

### Output

`seed_points.geojson` — a point layer in WGS84 with a `score` attribute (the original seed probability at each point location). Point density is controlled by `min_spacing_m`.

## Notes

- `createSeeds.R` and `RasterToVector.R` **must be run in the same R session**. `RasterToVector.R` uses the `seeds` `SpatRaster` object held in memory — it is not automatically re-loaded from disk.
- All intermediate rasters are removed via `rm()` + `gc()` in `createSeeds.R` to manage RAM on large DEMs.
- The Gaussian smoothing step uses `pad = TRUE` to avoid border artefacts; a subsequent `ifel` mask removes the zero-padded rim.
- Otsu thresholding operates on an 8-bit rescaled histogram (0–255) for algorithm consistency, then back-transforms to real-value units.