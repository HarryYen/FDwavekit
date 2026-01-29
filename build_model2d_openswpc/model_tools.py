from typing import Tuple

import numpy as np
from matplotlib.path import Path as MplPath
from scipy.ndimage import gaussian_filter


def grid_centers(
    x_min: float,
    x_max: float,
    z_min: float,
    z_max: float,
    grid_size: float,
) -> Tuple[np.ndarray, np.ndarray]:
    x_arr = np.arange(x_min, x_max, grid_size) + grid_size / 2
    z_arr = np.arange(z_min, z_max, grid_size) + grid_size / 2
    return x_arr, z_arr


def closest_index_in_array(matrix: np.ndarray, target_value: float) -> int:
    diff = np.abs(matrix - target_value)
    index = np.unravel_index(np.argmin(diff), diff.shape)
    return index[0]


def polygon_path(points: np.ndarray) -> MplPath:
    return MplPath(points)


def polygon_mask(path: MplPath, nx: int, nz: int) -> np.ndarray:
    mask = np.zeros((nx, nz), dtype=bool)
    for x_ii, z_ii in np.ndindex(nx, nz):
        if path.contains_point((x_ii, z_ii)):
            mask[x_ii, z_ii] = True
    return mask


def smooth_mask(mask: np.ndarray, sigma: float) -> np.ndarray:
    return gaussian_filter(mask.astype(float), sigma=sigma) > 0.5


def apply_mask(model: np.ndarray, mask: np.ndarray, new_value: float, use_percentage: bool) -> None:
    if use_percentage:
        model[mask] *= (1 + new_value / 100)
    else:
        model[mask] = new_value


def rectangle_path(
    x_min: float,
    z_min: float,
    x_max: float,
    z_max: float,
    grid_size: float,
    origin_x: float,
    origin_z: float,
) -> MplPath:
    x_lo, x_hi = sorted((x_min, x_max))
    z_lo, z_hi = sorted((z_min, z_max))
    x_list = np.array([x_lo, x_hi, x_hi, x_lo])
    z_list = np.array([z_lo, z_lo, z_hi, z_hi])
    x_ii_list = (x_list - origin_x) // grid_size + 1
    z_ii_list = (z_list - origin_z) // grid_size + 1
    x1_ii, x2_ii, x3_ii, x4_ii = x_ii_list
    z1_ii, z2_ii, z3_ii, z4_ii = z_ii_list
    return MplPath([(x1_ii, z1_ii), (x2_ii, z2_ii), (x3_ii, z3_ii), (x4_ii, z4_ii)])
