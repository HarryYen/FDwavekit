#%%
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import model_modules
import time

start = time.time()
# the basic model
nx, ny, nz = 10000, 34, 1050
initial_vel = 3 # km/s
grid_size  = 0.001 # km
x_min, y_min, z_min = 0, 0, -0.05

# automatically calculated
x_max = grid_size * nx + x_min
y_max = grid_size * ny + y_min
z_max = grid_size * nz + z_min

if __name__ == '__main__':

    # create the 2D model
    plane_model = np.ones((nx, nz)) * initial_vel

    ## create a fault zone
    plane_model = model_modules.set_up_fault_zone(fault_x   = 4.7, 
                                                  fault_dip = 60, 
                                                  fault_Vp  = 2.,
                                                  app_thick = 0.02,
                                                  model_arr = plane_model
    )
   
    ### create a rectangle area
    plane_model = model_modules.set_up_rectangle(
        arr       = plane_model, 
        x1 = 7.5, z1 = -0.05, x2 = 10, z2 = -0.05, x3 = 10, z3 = 0.2, x4 = 8, z4 = 0.2,  
        new_value = 1.5
    )
   
    ## padding 3 columns to the array for OpenSWPC reading
    plane_model = np.pad(plane_model, ((3,3),(3,3)), mode='edge')
    ## extending to 3D
    stacked_matrix = np.empty((nx+6, ny+6, nz+6))
    for y_ii in range(ny+6):
        print(y_ii)
        stacked_matrix[:,y_ii,:] = plane_model
    
    # model_modules.plot_2D_figures(stacked_matrix, 5., 0.125, 0)

    print('writing file...')
    model_modules.output(stacked_matrix)
 
    
    end = time.time()
    print(f'done! total time = {end-start:.1f}s')
    
# %%
