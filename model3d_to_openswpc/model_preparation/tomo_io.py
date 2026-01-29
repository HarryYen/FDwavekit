import numpy as np
import pandas as pd


def tomo_to_array(tomo_file, column_used):
    with open(tomo_file, "r") as f:
        line = f.readline()
        info = line.split()
        raw_nx = int(info[1])
        raw_ny = int(info[2])
        raw_nz = int(info[3])
    df = pd.read_csv(tomo_file, delimiter=r"\s+", skiprows=1, header=0)
    df_filter = df[["lon", "lat", "dep", column_used]]
    data_arr = np.zeros([raw_nx, raw_ny, raw_nz])

    lon_uniq = np.unique(np.array(df_filter["lon"]))
    lat_uniq = np.unique(np.array(df_filter["lat"]))
    dep_uniq = np.unique(np.array(df_filter["dep"]))

    index = 0
    for xx in range(raw_nx):
        for yy in range(raw_ny):
            for zz in range(raw_nz):
                data = df_filter.loc[index][column_used]
                data_arr[xx, yy, zz] = data
                index += 1
    return data_arr, lon_uniq, lat_uniq, dep_uniq
