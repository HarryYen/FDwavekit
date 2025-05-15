#%%
import xarray as xr
import pygmt
import sys

def extract_nc(nc_file):
    data = xr.open_dataset(nc_file)
    mu_arr = data['mu']
    lambda_arr = data['lambda']
    rho_arr = data['rho']
    return mu_arr, lambda_arr, rho_arr

def lame2vpvs(mu, lamb, rho):
    vp = ((lamb + 2*mu)/rho)**0.5
    vs = (mu/rho)**0.5
    return vp, vs

def pygmt_begin():
    fig = pygmt.Figure()
    return fig

def pygmt_plot(fig, arr):
    xmin, xmax = arr.x.min().item(), arr.x.max().item()
    zmin, zmax = arr.z.min().item(), arr.z.max().item()
    vmin, vmax = arr.min().item(), arr.max().item()
    fig.basemap(
        region=[xmin, xmax, zmin, zmax],
        projection='x1c/-1c',
        frame=["a1f1", "WSne", "x+lX(km)", "y+ldepth(km)"],
    )
    pygmt.makecpt(
        cmap = 'turbo',
        series = f'{vmin}/{vmax}',
        reverse = True,
    )
    fig.grdimage(
        grid = arr,
        cmap = True,
    )
    fig.colorbar(
        frame = ['x+l"Velocity (km/s)"'],
        position = 'JMR+o0.5c/0c+w3c/0.2c',
    )
    fig.show()

# ---------- PARAMETER ---------- #
nc_file = '/home/harry/Work/OpenSWPC-5.2.0/input_files_2d/MiDAS/ML_right_sf/out/ML_right_sf.xz.v.nc'
# --------------------------------#

if __name__ == '__main__':
    mu_arr, lambda_arr, rho_arr = extract_nc(nc_file)
    vp_arr, vs_arr = lame2vpvs(mu_arr, lambda_arr, rho_arr)
    print(vp_arr)
    print(vs_arr)
    for arr in [vp_arr, vs_arr]:
        fig = pygmt_begin()
        fig = pygmt_plot(fig, arr)
