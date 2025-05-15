#%%
from scipy.ndimage import gaussian_filter
from synmodel_creator import SynModelCreator

import numpy as np
import matplotlib.pyplot as plt
import sys
import time


def main():
    
    # ------------------------------------------------------------
    # Setting up the parameters 
    # (Please be the same as the one in the OpenSWPC)
    # ------------------------------------------------------------
    """
    Args:
        nx (int): number of grid points in x direction
        nz (int): number of grid points in z direction
        initial_vel (float): initial velocity in km/s
        grid_size (float): grid size in km
        x_min (float): minimum x coordinate in km
        z_min (float): minimum z coordinate in km
        specify_layered_model (bool): whether to specify a layered model (if False, the model will be a constant velocity model with initial_vel)
        layered_model_file (str): file name of the layered model (the format should be the same with ak135.csv)
    """
    config = {
        'nx': 2500,
        'nz': 3875,
        'initial_vel': 10, # km/s
        'grid_size': 0.04, # km
        'x_min': -40.0,
        'z_min': -5.0,
        'specify_layered_model': True,
        'layered_model_file': 'ak135.csv'
    }
    # ------------------------------------------------------------
    
    synmodel_creator = SynModelCreator(config)
    

    # --------------------------------------------------------- #
    # --------------- Initializing the array ------------------ #
    # --------------------------------------------------------- #
    synmodel_creator.create_homogeneous_model()
    
    # --------------------------------------------------------- #
    # ------------------ Design your model -------------------- #
    # --------------------------------------------------------- #
    if config['specify_layered_model']:
        # specify the layered model
        synmodel_creator.set_up_layered_structure()
    # --------------------------------------------------------- #
    # -------------specify your own structures----------------- #
    # --------------------------------------------------------- #
    synmodel_creator.set_up_fault_zone(surface_x = 5,
                                       dip_angle = 10, 
                                       deepest_z = 999., 
                                       app_thick = 999., 
                                       specified_value = 4., 
                                       percentage_flag = False)
    ### create a rectangle area
    # High-velocity boundary
    synmodel_creator.set_up_rectangle(
        x1 = 16, z1 = 20, x2 = 9999, z2 = 20, x3 = 9999, z3 = 9999, x4 = 16, z4 = 9999,  
        new_value = 7.2, percentage_flag = True
    )

    irregular_x = np.array([-15,  13., 12.5,    8,     5,     5,   -14, -16,  -14])
    irregular_z = np.array([42.5, 42.5, 60.5, 67.5,  75.5, 98.5, 98.5,  60, 42.5])
    irregular_z = irregular_z - 17.5
    points = list(zip(irregular_x, irregular_z))
    synmodel_creator.set_up_irregular_zone(points=points, 
                                           new_value=11.3, percentage_flag=True, 
                                           smooth_sigma=10, num_interpolated_points=100,
                                           s=1000, k=3)
   
   
    # Gaussian to blur the boundary
    synmodel_creator.blur_model(mode='nearest', sigma=50)
    # --------------------------------------------------------- #
    # -----------------------Post-Processing------------------- #
    # --------------------------------------------------------- #
    
    
    # padding 3 columns to the array for OpenSWPC reading (should NOT be ignored!)
    synmodel_creator.finalize_model()
    
    print('writing file...')
    synmodel_creator.output_2d()
    synmodel_creator.visualizing_array(vmin=1, vmax=10, label_name='Vp')

    
 


if __name__ == '__main__':
    start = time.time()
    # ------------------------- 
    main()
    # -------------------------
    end = time.time()
    print(f'done! total time = {end-start:.1f}s')
# %%
