# pip install geopandas rasterio shapely scikit-learn networkx pyproj pyogrio

import geopandas as gpd
import numpy as np
import rasterio
import networkx as nx
from shapely.geometry import LineString
from sklearn.neighbors import NearestNeighbors
from create_line_helper import sample_profile_ok

# -------------------------
# USER INPUTS
# -------------------------
pts_path = r"C:/Users/mimam/Dropbox/lider_data/output/points_created_v1.geojson"
dem_path = r"C:/Users/mimam/Dropbox/lider_data/input/pine_lake.tif"
out_lines_path = r"C:/Users/mimam/Dropbox/lider_data/output/connected_lines_v6.gpkg"     # or .geojson

# tuning
radius_m = 8      # search radius for neighbors
k_neighbors = 12     # max neighbors per point
alpha_dz = 10.0       # weight for elevation difference
min_len_m = 2.0     # remove tiny stubs
max_grade_edge = 0.29   # reject edges steeper than 25% grade (tune 0.10–0.35)
distance_weight = 1.0 # weight for horizontal distance in edge cost

# ---------- 1. Load points ----------
pts = gpd.read_file(pts_path)
if pts.crs is None:
    pts.set_crs(4326, inplace=True)
utm = pts.estimate_utm_crs()
pts_m = pts.to_crs(utm)

# ---------- 2. Sample DEM elevation ----------
with rasterio.open(dem_path) as src:
    pts_dem = pts_m.to_crs(src.crs)
    coords = [(g.x, g.y) for g in pts_dem.geometry]
    elev = np.array([v[0] for v in src.sample(coords)], dtype="float64")
pts_m["elev"] = elev
pts_m = pts_m[np.isfinite(pts_m["elev"])].copy()
pts_m.reset_index(drop=True, inplace=True)


# ---------- 3. Build weighted neighbor graph with in-between checks ----------
XY = np.array([(g.x, g.y) for g in pts_m.geometry])
nbrs = NearestNeighbors(n_neighbors=min(k_neighbors, len(pts_m)))
nbrs.fit(XY)
distances, indices = nbrs.kneighbors(XY)

G = nx.Graph(); G.add_nodes_from(range(len(pts_m)))

elev_arr = pts_m["elev"].to_numpy()
# Open DEM once more in this CRS for profile sampling
dem_src = rasterio.open(dem_path)  # your DEM (already used earlier)

for i in range(len(pts_m)):
    inds = np.asarray(indices[i]).ravel().astype(int)
    dsts = np.asarray(distances[i]).ravel()
    for j, d in zip(inds, dsts):
        if i >= j or d == 0:
            continue

        # coords in METERS (same CRS as DEM)
        p1 = (XY[i, 0], XY[i, 1])
        p2 = (XY[j, 0], XY[j, 1])

        # 1) reject if intermediate profile ever exceeds 12% grade
        if not sample_profile_ok(p1, p2, dem_src, max_grade=max_grade_edge, step_m=radius_m):
            continue

        # 3) still OK → compute simple weight
        dz = abs(float(elev_arr[i]) - float(elev_arr[j]))
        w  =distance_weight * float(d) + alpha_dz * dz
        # w  = alpha_dz * dz
        G.add_edge(i, j, weight=w, dist=float(d), dz=float(dz))


# ---------- 4. Minimum Spanning Forest ----------
H = nx.Graph()
for comp in nx.connected_components(G):
    sub = G.subgraph(comp)
    mst = nx.minimum_spanning_tree(sub, weight="weight")
    H = nx.compose(H, mst)

# ---------- 5. Decompose forest into paths ----------
def decompose_into_paths(graph: nx.Graph):
    paths = []
    used = set()
    for n in graph.nodes():
        if graph.degree(n) != 2:
            for nb in graph.neighbors(n):
                e = (min(n, nb), max(n, nb))
                if e in used:
                    continue
                path = [n, nb]
                used.add(e)
                prev, cur = n, nb
                while graph.degree(cur) == 2:
                    nbs = list(graph.neighbors(cur))
                    nxt = nbs[0] if nbs[0] != prev else nbs[1]
                    e2 = (min(cur, nxt), max(cur, nxt))
                    if e2 in used:
                        break
                    path.append(nxt)
                    used.add(e2)
                    prev, cur = cur, nxt
                paths.append(path)
    # cycles
    for u, v in graph.edges():
        e = (min(u, v), max(u, v))
        if e in used:
            continue
        path = [u, v]
        used.add(e)
        prev, cur = u, v
        while True:
            nbs = list(graph.neighbors(cur))
            nxt = nbs[0] if nbs[0] != prev else nbs[1]
            e2 = (min(cur, nxt), max(cur, nxt))
            if e2 in used:
                break
            path.append(nxt)
            used.add(e2)
            prev, cur = cur, nxt
        paths.append(path)
    return paths

paths = decompose_into_paths(H)

# ---------- 6. Convert to LineStrings ----------
def paths_to_lines(paths, gdf_points_metric_crs):
    lines = []
    attrs = []
    XY = np.array([(geom.x, geom.y) for geom in gdf_points_metric_crs.geometry])
    for node_ids in paths:
        if len(node_ids) < 2:
            continue
        coords = [(XY[i, 0], XY[i, 1]) for i in node_ids]
        lines.append(LineString(coords))
        attrs.append({"n_points": len(node_ids)})
    return gpd.GeoDataFrame(attrs, geometry=lines, crs=gdf_points_metric_crs.crs)

lines_gdf = paths_to_lines(paths, pts_m)
lines_gdf = lines_gdf[lines_gdf.length >= min_len_m].copy()

# ---------- 7. Save ----------
lines_gdf.to_file(out_lines_path, driver="GPKG")
print(f"Saved {len(lines_gdf)} lines → {out_lines_path}")