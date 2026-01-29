import numpy as np


def len_to_grid_number(beg, end, dx):
    span = (end - beg)
    nx = span // dx
    return int(nx)


def arange_z_arr(nz, dz, z_beg):
    z_beg_new = round(z_beg + dz / 2 - dz * 3, 3)
    z_max = round(z_beg + nz * dz - dz / 2 + dz * 3, 3)
    z_arr = np.arange(z_beg_new, z_max + dz / 2, dz)
    print(z_beg_new, z_max)
    return z_arr
