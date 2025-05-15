import numpy as np
import pandas as pd
import pygmt
import itertools
import multiprocessing
import time
import yaml
from scipy.interpolate import RegularGridInterpolator
from functools import partial

def load_config(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def len_to_grid_number(beg, end, dx):
    len = (end - beg)
    nx = len//dx
    return int(nx)

def tomo_to_array(tomo_file, column_used):
    with open(tomo_file, 'r') as f:
        line = f.readline()
        info = line.split()
        raw_nx = int(info[1])
        raw_ny = int(info[2])
        raw_nz = int(info[3])
    df = pd.read_csv(tomo_file, delimiter=r'\s+', skiprows=1, header=0)
    df_filter = df[['lon','lat','dep',column_used]]
    data_arr = np.zeros([raw_nx,raw_ny,raw_nz])

    lon_uniq = np.unique(np.array(df_filter['lon']))
    lat_uniq = np.unique(np.array(df_filter['lat']))
    dep_uniq = np.unique(np.array(df_filter['dep']))

    index = 0
    for xx in range(raw_nx):
        for yy in range(raw_ny):
            for zz in range(raw_nz): 
                data = df_filter.loc[index][column_used]
                data_arr[xx,yy,zz] = data
                index += 1
    return data_arr, lon_uniq, lat_uniq, dep_uniq
    

def project_hor_and_vert_lines(cx, cy, nx, ny, dx, dy, x_half_len, y_half_len):
    x_starting_len = -(x_half_len) + dx/2 - 3*dx
    x_ending_len = -(x_half_len) + dx*nx - dx/2 + 3*dx
    y_starting_len = -(y_half_len) + dy/2 - 3*dy
    y_ending_len = -(y_half_len) + dy*ny - dy/2 + 3*dy

    df_x = pygmt.project(center=[cx, cy], azimuth=90, length=[x_starting_len, x_ending_len],
                                 generate=dx, unit=True)
    df_y = pygmt.project(center=[cx, cy], azimuth=0, length=[y_starting_len, y_ending_len],
                                 generate=dy, unit=True)
    x_proj_arr = np.array(df_x.r)
    y_proj_arr = np.array(df_y.s)

    return x_proj_arr, y_proj_arr

def arange_z_arr(nz, dz, z_beg):
    z_beg_new = round(z_beg + dz/2 - dz*3, 3)
    z_max = round(z_beg + nz*dz - dz/2 + dz*3, 3)
    z_arr = np.arange(z_beg_new, z_max+dz/2, dz) # z_max+dz for avoiding step stop before z_max
    print(z_beg_new, z_max)
    return z_arr

 

def parallel_create_final_arr(interp_func, x_arr, y_arr, z_arr, params):
    ii, jj, kk = params[0], params[1], params[2]
    value = interp_func([x_arr[ii], y_arr[jj], z_arr[kk]])[0]

    return value
    

def output(arr, specify_data):
    f = open(f'3d_model_{specify_data}.dat', 'w+')
    for ii in range(len(arr)):
        f.write(f'{arr[ii]:.3f}\n')
    f.close()
