# pip install rasterio geopandas shapely scikit-image scikit-learn numpy pyogrio

import numpy as np
import rasterio
import geopandas as gpd
from rasterio.transform import xy
from skimage.morphology import (
    binary_opening, skeletonize, remove_small_objects, footprint_rectangle
)
from skimage.measure import label, regionprops
from sklearn.neighbors import NearestNeighbors

# ------------------ USER INPUTS ------------------
prob_tif   = r"C:/Users/mimam/Dropbox/lider_data/output/seeds_v2.tif"  # continuous road-likelihood raster (0..1 preferred)
out_points = r"C:/Users/mimam/Dropbox/lider_data/output/seed_points_monotonic.gpkg"

# Threshold mode (choose ONE)
prob_threshold  = 0.25      # absolute threshold

# Light cleanup on the binary mask (before skeleton)
do_opening = False          # opening may break thin roads; start False
opening_size = 3            # if do_opening=True, 3x3 rectangle

# Skeleton post-clean
min_skel_blob_px = 1       # drop tiny skeleton fragments (in pixels)

# Keep long corridors only (monotonic wrt threshold)
min_component_len_m = 1.0  # minimum skeleton component length to keep (meters)

# Point thinning (meters) for even spacing along roads
min_spacing_m = 8.0

# ------------------ LOAD & NORMALIZE ------------------
with rasterio.open(prob_tif) as src:
    prob = src.read(1).astype("float32")
    transform, crs = src.transform, src.crs
    nodata = src.nodata
    # pixel size in meters (assumes near-square pixels)
    pix_m = float(np.hypot(src.transform.a, src.transform.e))

if nodata is not None:
    prob = np.where(prob == nodata, np.nan, prob)

finite = np.isfinite(prob)
if not finite.any():
    raise RuntimeError("No finite pixels in seeds raster.")

# ------------------ THRESHOLD (MONOTONIC) ------------------
thr = float(prob_threshold)

mask = (prob >= thr) & finite
print(f"[A] threshold={thr:.3f}  → mask pixels: {int(mask.sum())}")

# optional: very light opening on the binary mask
if do_opening:
    mask = binary_opening(mask, footprint_rectangle(shape=(opening_size, opening_size)))
    print(f"[A] after opening: {int(mask.sum())} pixels")

# ------------------ SKELETONIZE FIRST ------------------
skel = skeletonize(mask)
print(f"[B] skeleton pixels: {int(skel.sum())}")

# remove tiny skeleton fragments (in pixels)
if min_skel_blob_px > 0:
    skel = remove_small_objects(skel, min_size=min_skel_blob_px)
    print(f"[B] skeleton after small-object removal: {int(skel.sum())} px")

# ------------------ KEEP LONG COMPONENTS ------------------
lab = label(skel, connectivity=2)
props = regionprops(lab)

keep_labels = {rp.label for rp in props if (rp.area * pix_m) >= min_component_len_m}
if not keep_labels:
    # fallback: keep top few longest components so we don't return empty
    top = sorted(props, key=lambda r: r.area, reverse=True)[:5]
    keep_labels = {r.label for r in top}
    print("[C] no component ≥ min length; keeping top 5 longest as fallback.")

skel_keep = np.isin(lab, list(keep_labels))
print(f"[C] kept skeleton pixels (after length filter): {int(skel_keep.sum())}")

# ------------------ RASTER → POINTS ------------------
rows, cols = np.where(skel_keep)
if len(rows) == 0:
    raise RuntimeError("No skeleton pixels survived. Lower threshold or relax filters.")

xs, ys = rasterio.transform.xy(transform, rows, cols)
scores = prob[rows, cols].astype("float32")

gdf = gpd.GeoDataFrame({"score": scores},
                       geometry=gpd.points_from_xy(xs, ys),
                       crs=crs)

# ------------------ THIN POINTS (EVEN SPACING IN METERS) ------------------
utm = gdf.estimate_utm_crs()
gdf_m = gdf.to_crs(utm)
XY = np.c_[gdf_m.geometry.x, gdf_m.geometry.y]

# Greedy thinning: keep highest score first, suppress neighbors within min_spacing_m
order = np.argsort(-gdf_m["score"].to_numpy())
alive = np.ones(len(gdf_m), dtype=bool)
keep_idx = []
nbr = NearestNeighbors(radius=min_spacing_m).fit(XY)

for i in order:
    if not alive[i]:
        continue
    keep_idx.append(i)
    idx = nbr.radius_neighbors([XY[i]], return_distance=False)[0]
    alive[idx] = False
    alive[i] = True

gdf_thin = gdf_m.iloc[keep_idx].to_crs(crs).copy()
print(f"[D] thinned points: {len(gdf_thin)}")

# ------------------ SAVE ------------------
gdf_thin.to_file(out_points, driver="GPKG")
print(f"Saved {len(gdf_thin)} points → {out_points}")
