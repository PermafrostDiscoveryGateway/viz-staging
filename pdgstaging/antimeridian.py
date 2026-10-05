"""Repair polygons that cross the antimeridian after reprojection to lon/lat.

Reprojecting a polygon that straddles +/-180 degrees (for example from a polar
projection) gives a ring whose vertices jump between about +180 and -180. In
plain lon/lat such a ring spans the whole globe: it is usually invalid, its
area and centroid are wrong, and it is drawn as a band around the world.

The fix is to unwrap the ring (add 360 to its negative longitudes so it is
contiguous), take area and centroid from that contiguous shape, and then split
it at 180 degrees so the stored geometry stays inside [-180, 180].
"""

from __future__ import annotations

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon, box

# Longitude jump between consecutive vertices that marks a wrapped ring.
JUMP_DEGREES = 180.0


def _rings(geometry):
    for part in shapely.get_parts(geometry):
        if isinstance(part, Polygon) and not part.is_empty:
            yield part.exterior
            yield from part.interiors


def crosses_antimeridian(geometry):
    """Return True if any ring has consecutive vertices more than 180 degrees apart in longitude."""
    for ring in _rings(geometry):
        x = np.asarray(ring.coords)[:, 0]
        if x.size > 1 and np.abs(np.diff(x)).max() > JUMP_DEGREES:
            return True
    return False


def unwrap_antimeridian(geometry):
    """Return ``geometry`` with negative longitudes shifted by +360 if it crosses.

    Geometry that does not cross the antimeridian, and geometry that is still
    invalid after unwrapping (for example a polygon enclosing a pole), is
    returned unchanged (the same object).
    """
    if geometry is None or geometry.is_empty or not crosses_antimeridian(geometry):
        return geometry

    def shift(coords):
        coords = np.array(coords, dtype=float)
        coords[coords[:, 0] < 0, 0] += 360.0
        return coords

    unwrapped = shapely.transform(geometry, shift)
    return unwrapped if unwrapped.is_valid else geometry


def split_at_antimeridian(geometry):
    """Cut an unwrapped geometry at 180 degrees and fold the east part back to [-180, 180].

    Geometry that does not extend past 180 degrees is returned unchanged.
    """
    if geometry.bounds[2] <= 180.0:
        return geometry

    west_of_180 = shapely.intersection(geometry, box(-180.0, -90.0, 180.0, 90.0))
    past_180 = shapely.intersection(geometry, box(180.0, -90.0, 540.0, 90.0))
    folded = shapely.transform(past_180, lambda c: c - [360.0, 0.0])

    parts = [
        part
        for piece in (west_of_180, folded)
        for part in shapely.get_parts(piece)
        if isinstance(part, Polygon) and not part.is_empty and part.area > 0
    ]
    return MultiPolygon(parts)


def wrap_longitude(x):
    """Wrap longitudes into [-180, 180)."""
    return (np.asarray(x, dtype=float) + 180.0) % 360.0 - 180.0
