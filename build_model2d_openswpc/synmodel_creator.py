from pathlib import Path
from typing import Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MplPath
from scipy.interpolate import splprep, splev
from scipy.ndimage import gaussian_filter


class SynModelCreator:
    """Build and edit a 2D velocity model for OpenSWPC."""

    def __init__(self, config: dict) -> None:
        self.nx = config['nx']
        self.nz = config['nz']
        self.initial_vel = config['initial_vel']
        self.grid_size = config['grid_size']
        self.x_min = config['x_min']
        self.z_min = config['z_min']
        self.specify_layered_model = config['specify_layered_model']
        self.input_dir = Path(config.get('input_dir', '.'))
        self.output_dir = Path(config.get('output_dir', '.'))
        self.output_dir.mkdir(parents=True, exist_ok=True)

        layered_model_file = config['layered_model_file']
        layered_model_path = Path(layered_model_file)
        if not layered_model_path.is_absolute():
            candidate_path = self.input_dir / layered_model_file
            if candidate_path.exists():
                layered_model_path = candidate_path
        self.layered_model_file = str(layered_model_path)

        # init 
        self.get_xz_max()
        self.plane_model = np.ones((self.nx, self.nz))
    
    def get_xz_max(self) -> None:
        self.x_max = self.grid_size * self.nx + self.x_min
        self.z_max = self.grid_size * self.nz + self.z_min

    def _grid_centers(self) -> Tuple[np.ndarray, np.ndarray]:
        x_arr = np.arange(self.x_min, self.x_max, self.grid_size) + self.grid_size / 2
        z_arr = np.arange(self.z_min, self.z_max, self.grid_size) + self.grid_size / 2
        return x_arr, z_arr
    
    def create_homogeneous_model(self) -> None:
        """
        Create an initial homogeneous model with the given parameters.
        """
        self.plane_model *= self.initial_vel

    def set_up_layered_structure(self) -> None:
        """
        Set up the initial model based on the specified 1D layer model.
        """
        vel_data = np.loadtxt(self.layered_model_file)
        dep_values = vel_data[:, 0]
        vel_values = vel_data[:, 1]
        nonzero_min = np.min(vel_values[np.nonzero(vel_values)])
        
        for ii, dep in enumerate(dep_values):
            if dep >= self.z_max:
                break
            vel = vel_values[ii]
            if vel == 0.0:
                vel = nonzero_min
            dep_ii = int(np.ceil(round(dep - self.z_min, 4) / self.grid_size))
            dep_ii = max(dep_ii, 0)
            self.plane_model[:, dep_ii:] = vel

    def set_up_rectangle(
        self,
        x_min: float,
        z_min: float,
        x_max: float,
        z_max: float,
        new_value: float,
        use_percentage: bool,
    ) -> None:
        """
        Add a rectangle area in the model.
        Args:
            x_min, z_min, x_max, z_max (float): coordinates of two opposite corners
            new_value (float): value to set in the rectangle area
            use_percentage (bool): if True, set the specifed value as a percentage of the original velocity
                                    True: NEW = OLD * (1 + new_value/100)
                                    False: NEW = new_value
        """
        x_lo, x_hi = sorted((x_min, x_max))
        z_lo, z_hi = sorted((z_min, z_max))
        x_list = np.array([x_lo, x_hi, x_hi, x_lo])
        z_list = np.array([z_lo, z_lo, z_hi, z_hi])
        x_ii_list = (x_list - self.x_min) // self.grid_size + 1
        z_ii_list = (z_list - self.z_min) // self.grid_size + 1
        x1_ii, x2_ii, x3_ii, x4_ii = x_ii_list
        z1_ii, z2_ii, z3_ii, z4_ii = z_ii_list
        p = MplPath([(x1_ii, z1_ii), (x2_ii, z2_ii), (x3_ii, z3_ii), (x4_ii, z4_ii)])
        scale = 1 + new_value / 100
        for x_ii in range(self.nx):
            for z_ii in range(self.nz):
                is_in_rectangle = p.contains_point((x_ii, z_ii))
                if is_in_rectangle:

                    if use_percentage:
                        self.plane_model[x_ii, z_ii] = self.plane_model[x_ii, z_ii] * scale
                    else:
                        self.plane_model[x_ii, z_ii] = new_value

    @staticmethod
    def find_closest_index_in_array(matrix: np.ndarray, target_value: float) -> int:
        diff = np.abs(matrix - target_value)
        index = np.unravel_index(np.argmin(diff), diff.shape)
        return index[0]

    def set_up_fault_zone(
        self,
        surface_x: float,
        dip_angle: float,
        deepest_z: float,
        app_thick: float,
        specified_value: float,
        use_percentage: bool,
    ) -> None:
        """
        Set up a fault zone with the given parameters.
        Args:
            surface_x (float): x coordinate of the intersection between free surface and fault line
            dip_angle (float): dip angle of the fault in degrees
            deepest_z (float): deepest z coordinate of the fault
            app_thick (float): apparent thickness of the fault zone
            specified_value (float): value to set in the fault zone
            use_percentage (bool): if True, set the specifed value as a percentage of the original velocity
                                    True: NEW = OLD * (1 + specified_value/100)
                                    False: NEW = specified_value
        """
        slope = np.tan(np.deg2rad(dip_angle))

        x_arr, z_arr = self._grid_centers()
        z_arr_index = range(self.nz)
        fault_top = slope * (-x_arr + surface_x)
        fault_bot = fault_top + app_thick
        root1_arr = -1 * fault_top  
        root2_arr = -1 * fault_bot  

        for ii in range(len(x_arr)):
            root_min = min(root1_arr[ii], root2_arr[ii])
            root_max = max(root1_arr[ii], root2_arr[ii])
            if root_max > deepest_z:
                if root_min > deepest_z: 
                    continue
                elif root_min <= deepest_z:
                    root_max = deepest_z
            
            if not ((root_min <= self.z_min and root_max <= self.z_min) or (root_min >= self.z_max and root_max >= self.z_max)):
                root_min_ii = self.find_closest_index_in_array(z_arr, root_min)
                root_max_ii = self.find_closest_index_in_array(z_arr, root_max)
                
                if (root_min_ii == z_arr_index[0] and root_max_ii == z_arr_index[0]) or (root_min_ii == z_arr_index[-1] and root_max_ii == z_arr_index[-1]):
                    continue
                else:
                    if use_percentage:
                        self.plane_model[ii, root_min_ii:root_max_ii] = self.plane_model[ii, root_min_ii:root_max_ii] * (1 + specified_value / 100)
                    else:    
                        self.plane_model[ii, root_min_ii:root_max_ii] = specified_value
    

    def set_up_irregular_zone(
        self,
        points: Sequence[Tuple[float, float]],
        new_value: float,
        use_percentage: bool,
        s: float,
        k: int,
        num_interpolated_points: int = 100,
    ) -> None:
        """
        Add an irregular polygonal zone and assign a new value inside it.

        The boundary is created from control points using a closed B-spline,
        then rasterized into a mask and optionally smoothed.

        Args:
            points (Sequence[Tuple[float, float]]): Control points in (x, z) coordinates.
            new_value (float): Fixed value to set, or percentage change.
            use_percentage (bool): If True, apply NEW = OLD * (1 + new_value / 100).
                                   If False, apply NEW = new_value.
            s (float): Spline smoothness for `splprep` (0 = exact fit, larger = smoother).
            k (int): Spline degree for `splprep` (typically 1-5).
            num_interpolated_points (int): Number of points used to densify the spline.
        """
        
        # Convert (x, z) coordinates to index coordinates
        points_arr = np.asarray(points)
        x_indices = ((points_arr[:, 0] - self.x_min) // self.grid_size).astype(int)
        z_indices = ((points_arr[:, 1] - self.z_min) // self.grid_size).astype(int)
        
        # Use B-spline interpolation to smooth the boundary with adjustable smoothness
        tck, u = splprep([x_indices, z_indices], s=s, k=k, per=True)
        u_fine = np.linspace(0, 1, num_interpolated_points)
        x_smooth, z_smooth = splev(u_fine, tck)
        
        # Create region mask
        polygon_path = MplPath(np.column_stack([x_smooth, z_smooth]))
        mask = np.zeros_like(self.plane_model, dtype=bool)
        
        for x_ii, z_ii in np.ndindex(self.nx, self.nz):
            if polygon_path.contains_point((x_ii, z_ii)):
                mask[x_ii, z_ii] = True
        
        # Apply new values
        if use_percentage:
            self.plane_model[mask] *= (1 + new_value / 100)
        else:
            self.plane_model[mask] = new_value
    
    def blur_model(self, mode: str = 'nearest', sigma: float = 5) -> None:
        """
        Apply Gaussian smoothing to the model.
        
        Args:
            mode (str): Mode for Gaussian filter. Default is 'nearest'.
            sigma (float): Standard deviation for Gaussian kernel.
        """
        self.plane_model = gaussian_filter(self.plane_model, mode=mode, sigma=sigma)

    def finalize_model(self) -> None:
        """
        Now we need to pad the model with 3 grids each side because of the setting in OpenSWPC
        """
        self.plane_model = np.pad(self.plane_model, ((3, 3), (3, 3)), mode='edge')
        self.nx += 6
        self.nz += 6

    def output_2d(self) -> None:
        output_path = self.output_dir / 'model_Vp_2d.dat'
        print(self.nx, self.nz)
        with open(output_path, 'w+') as f:
            for ii in range(self.nx):
                for kk in range(self.nz):
                    f.write(f'{self.plane_model[ii, kk]:.3f}\n')

    def visualizing_array(self, vmin: float, vmax: float, label_name: str) -> None:
        fig, _ = plt.subplots(1, 1, figsize=(5, 5))
        im = plt.imshow(
            self.plane_model[:, :].T,
            cmap='coolwarm_r',
            extent=[self.x_min, self.x_max, self.z_max, self.z_min],
            vmin=vmin,
            vmax=vmax,
        )
        plt.xlabel('X(km)')
        plt.ylabel('Z(km)')
        cbar = plt.colorbar(im, shrink=0.4)
        cbar.set_label(label_name)
        fig.savefig(self.output_dir / 'XZ.png', dpi=300)
