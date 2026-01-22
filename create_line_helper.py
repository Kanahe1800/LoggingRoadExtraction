import math
from shapely.geometry import LineString
import numpy as np
from pyproj import Transformer

def sample_profile_ok(p1, p2, dem_src, max_grade=0.12, step_m=5.0):
    """
    Check the straight-line profile between p1=(x,y) and p2=(x,y) in DEM CRS.
    Returns True if all intermediate segment grades <= max_grade.
    """
    x1, y1 = p1; x2, y2 = p2
    dx, dy = x2 - x1, y2 - y1
    dist = math.hypot(dx, dy)
    if dist == 0:
        return False
    n_steps = max(1, int(dist // step_m))
    # coordinates to sample (including endpoints)
    xs = np.linspace(x1, x2, n_steps + 1)
    ys = np.linspace(y1, y2, n_steps + 1)
    elevs = [v[0] for v in dem_src.sample(list(zip(xs, ys)))]
    # if any nodata → reject
    if any((e is None) or (np.isnan(e)) for e in elevs):
        return False
    # check grade between consecutive samples
    for k in range(len(elevs) - 1):
        run = step_m if n_steps > 1 else dist
        rise = abs(float(elevs[k+1]) - float(elevs[k]))
        grade = rise / max(run, 1e-6)
        if grade > max_grade:
            return False
    return True
