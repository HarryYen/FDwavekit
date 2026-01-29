#%%
from dataclasses import asdict, dataclass
import time
import numpy as np
from synmodel_creator import SynModelCreator


@dataclass
class ModelConfig:
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
        input_dir (str): folder for initial 1D velocity model files
        output_dir (str): folder for model outputs
    """
    nx: int = 250
    nz: int = 388
    grid_size: float = 0.4  # km
    x_min: float = -40.0
    z_min: float = -5.0
    
    initial_vel: float = 10.0  # km/s
    specify_layered_model: bool = True
    layered_model_file: str = 'ak135.csv'
    
    input_dir: str = 'input'
    output_dir: str = 'output'


def main() -> None:
    
    # Model parameters (match your OpenSWPC settings).
    config = asdict(ModelConfig())
    
    synmodel_creator = SynModelCreator(config)

    # Initialize the array.
    synmodel_creator.create_homogeneous_model()

    # Design your model.
    if config['specify_layered_model']:
        synmodel_creator.set_up_layered_structure()

    # Specify your own structures.
    synmodel_creator.set_up_fault_zone(surface_x = -20,
                                       dip_angle = 20, 
                                       deepest_z = 40., 
                                       app_thick = 10., 
                                       specified_value = 8., 
                                       use_percentage = False)
    
    ### create a rectangle area
    synmodel_creator.set_up_rectangle(
        x_min = -20, z_min = 20, x_max = 20, z_max = 40,
        new_value = 7.2, use_percentage = False
    )

    # Demo: irregular zone (comment out if not needed).
    irregular_x = np.array([-20, -20, 20, 20])
    irregular_z = np.array([60, 80, 80, 60])
    points = list(zip(irregular_x, irregular_z))
    synmodel_creator.set_up_irregular_zone(
        points=points,
        new_value=-50,
        use_percentage=True,
        num_interpolated_points=100,
        s=3,
        k=3,
    )
   
   
    # Gaussian to blur the boundary.
    # synmodel_creator.blur_model(mode='nearest', sigma=5)

    # -----------------------------------------------------------------------------
    # Post-processing.
    # padding 3 columns to the array for OpenSWPC reading (should NOT be ignored!)
    # -----------------------------------------------------------------------------
    synmodel_creator.finalize_model()
    
    print('writing file...')
    synmodel_creator.output_2d()
    synmodel_creator.visualizing_array(vmin=1, vmax=10, label_name='Vp')

    
 


if __name__ == '__main__':
    start = time.perf_counter()
    # ------------------------- 
    main()
    # -------------------------
    elapsed = time.perf_counter() - start
    print(f'done! total time = {elapsed:.1f}s')
# %%
