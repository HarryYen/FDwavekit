#%%
import xarray as xr
import sys
import numpy as np
import pygmt
import os

# ----------------------------------------
# Parameters
# ----------------------------------------
work_dir = '../../OpenSWPC-5.2.0/input_files'
project_name = 'test2d_eq'
wavefield_type = 'Vz'
nc_f = os.path.join(work_dir, project_name, 'out', f'{project_name}.xz.v.nc')
mul = 5
# ----------------------------------------

    
saving_file = project_name 
region = [0, 10, -1, 2]    
ds = xr.open_dataset(nc_f)

vmin, vmax = ds[wavefield_type].actual_range
xmin, xmax = ds['x'].actual_range
zmin, zmax = ds['z'].actual_range

vmax_abs = np.max([np.abs(vmin), np.abs(vmax)])

pygmt.makecpt(cmap='polar', series=[-vmax_abs, vmax_abs], background=True)
grouped_T = ds.groupby('t')
Vs_model_grid = np.sqrt(ds['mu'] / ds['rho'])



region = [xmin, xmax, zmin, zmax]
region_show = region
print(region_show)

t_text_loc = [xmin + (xmax-xmin)*0.02, zmin + (zmax-zmin)*0.02]

for it, t_value in enumerate(ds.t.values):

    fig = pygmt.Figure()
    fig.basemap(
        region     = region,
        projection = 'x0.25c/-0.25c',
        frame=["WSne", "xa10f10", "ya10f10","x+lX(km)","y+lDepth(km)"],
    )
    print(f't={t_value}')
    
    
    vz_data = ds[wavefield_type].isel(t=it).load()
    vz_data = vz_data * mul
    print(vz_data)
    

    fig.grdimage(
        grid = vz_data,
        cmap = True,
        transparency = 0,
        nan_transparent = True,
        frame = ["news"],
        region = region_show,
        
    )

    fig.grdcontour(grid=Vs_model_grid, interval=0.05, 
                transparency=90, pen='0.01p,black',
    )


    fig.text(text = f'T={t_value:.3f}s',
            x = t_text_loc[0],
            y = t_text_loc[1],
            fill = 'gray',
            justify = 'CM',
            font = '14p,Helvetica-Bold'
    )
    
    
    os.makedirs(f'snapshots/{saving_file}', exist_ok=True)
    fig.savefig(f'snapshots/{saving_file}/{t_value:.5f}.png', dpi=300)

            
# %%
