import numpy as np
import pygmt


def project_hor_and_vert_lines(cx, cy, nx, ny, dx, dy, x_half_len, y_half_len):
    x_starting_len = -(x_half_len) + dx / 2 - 3 * dx
    x_ending_len = -(x_half_len) + dx * nx - dx / 2 + 3 * dx
    y_starting_len = -(y_half_len) + dy / 2 - 3 * dy
    y_ending_len = -(y_half_len) + dy * ny - dy / 2 + 3 * dy

    df_x = pygmt.project(
        center=[cx, cy], azimuth=90, length=[x_starting_len, x_ending_len], generate=dx, unit=True
    )
    df_y = pygmt.project(
        center=[cx, cy], azimuth=0, length=[y_starting_len, y_ending_len], generate=dy, unit=True
    )
    x_proj_arr = np.array(df_x.r)
    y_proj_arr = np.array(df_y.s)

    return x_proj_arr, y_proj_arr
