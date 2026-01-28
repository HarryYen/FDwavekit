#%%
import pygmt
import itertools
import multiprocessing
import time
import numpy as np
import pandas as pd
import sys

from functools import partial
from scipy.interpolate import RegularGridInterpolator
from model_preparation.module_interp3d_2 import *


if __name__ == '__main__':
    # ----------LOAD PARAMETERS FROM INPUT FILE----------
    config = load_config('input.yaml')
    tomo_file = config['tomo_file']
    dx = config['dx']
    dy = config['dy']
    dz = config['dz']
    cx = config['cx']
    cy = config['cy']
    x_half_len = config['x_half_len']
    y_half_len = config['y_half_len']
    z_beg = config['z_beg']
    z_end = config['z_end']
    interpolation_method = config['interpolation_method']
    vector = config['vector']
    default_value = config['default_value']
    cores_num = config['cores_num']
    # ---------------------------------------------------
    print(tomo_file) 
    nx = len_to_grid_number(beg=-x_half_len, end=x_half_len, dx=dx)
    ny = len_to_grid_number(beg=-y_half_len, end=y_half_len, dx=dy)
    nz = len_to_grid_number(beg=z_beg, end=z_end, dx=dz)
    tomo_dir = f'model_input/{tomo_file}'
    data_arr, lon_uniq, lat_uniq, dep_uniq = tomo_to_array(tomo_dir, vector)
    x_arr, y_arr = project_hor_and_vert_lines(cx=cx, cy=cy, nx=nx,ny=ny, dx=dx, dy=dy,
                                            x_half_len=x_half_len,y_half_len=y_half_len)
    z_arr = arange_z_arr(nz=nz, dz=dz, z_beg=z_beg)

    interp_nx, interp_ny, interp_nz = len(x_arr), len(y_arr), len(z_arr)
    print(interp_nx, interp_ny, interp_nz)
    paramlist = list(itertools.product(range(interp_nx),range(interp_ny),range(interp_nz)))


    specify_data = vector
    interpolating_function = RegularGridInterpolator(points=(lon_uniq, lat_uniq, dep_uniq), values=data_arr, 
                                                    method=interpolation_method, fill_value=default_value, bounds_error=False)

    #Generate processes equal to the number of cores
    with multiprocessing.Pool(cores_num) as pool:
        res = pool.map(partial(parallel_create_final_arr, interpolating_function, x_arr, y_arr, z_arr), paramlist)
    output(arr=res, specify_data=specify_data)
    print(f'{specify_data} DONE!!')
