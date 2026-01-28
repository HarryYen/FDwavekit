from matplotlib.path import Path
from scipy.interpolate import splprep, splev
from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt


class SynModelCreator:
    def __init__(self, config):
        self.nx = config['nx']
        self.nz = config['nz']
        self.initial_vel = config['initial_vel']
        self.grid_size = config['grid_size']
        self.x_min = config['x_min']
        self.z_min = config['z_min']
        self.specify_layered_model = config['specify_layered_model']
        self.layered_model_file = config['layered_model_file']

        # init 
        self.get_xz_max()
        self.plane_model = np.ones((self.nx, self.nz))
    
    def get_xz_max(self):
    
        self.x_max = self.grid_size * self.nx + self.x_min
        self.z_max = self.grid_size * self.nz + self.z_min
    
    def create_homogeneous_model(self):
        """
        Create an initial homogeneous model with the given parameters.
        """
        self.plane_model = self.plane_model * self.initial_vel

    def set_up_layered_structure(self):
        """
        Set up the initial model based on the specified 1D layer model.
        """
        vel_data = np.loadtxt(self.layered_model_file)
        dep_values = vel_data[:,0]
        vel_values = vel_data[:,1]
        nonzero_min = np.min(vel_values[np.nonzero(vel_values)])
        
        for ii, dep in enumerate(dep_values):
            if dep >= self.z_max:
                break
            vel = vel_values[ii]
            if vel == 0.0:
                vel = nonzero_min
            dep_ii = int(np.ceil(round(dep-self.z_min, 4) / self.grid_size))
            if dep_ii < 0:
                dep_ii = 0
            self.plane_model[:, dep_ii:] = vel

    def set_up_rectangle(self, x1, z1, x2, z2, x3, z3, x4, z4, new_value, percentage_flag):
        """
        add a rectange area in the model
        Args:
            x1, z1, x2, z2, x3, z3, x4, z4 (float): coordinates of the rectangle corners
            new_value (float): value to set in the rectangle area
            percentage_flag (bool): if True, set the specifed value as a percentage of the original velocity
                                    True: NEW = OLD * (1 + new_value/100)
                                    False: NEW = new_value
        """
        x_list = np.array([x1, x2, x3, x4])
        z_list = np.array([z1, z2, z3, z4])
        x_ii_list = (x_list - self.x_min) // self.grid_size + 1
        z_ii_list = (z_list - self.z_min) // self.grid_size + 1
        x1_ii, x2_ii, x3_ii, x4_ii = x_ii_list[0], x_ii_list[1], x_ii_list[2], x_ii_list[3]
        z1_ii, z2_ii, z3_ii, z4_ii = z_ii_list[0], z_ii_list[1], z_ii_list[2], z_ii_list[3]
        p = Path([(x1_ii, z1_ii), (x2_ii, z2_ii), (x3_ii, z3_ii), (x4_ii, z4_ii)])
        for x_ii in range(self.nx):
            for z_ii in range(self.nz):
                is_in_rectangle = p.contains_point((x_ii, z_ii))
                if is_in_rectangle:

                    if percentage_flag:
                        self.plane_model[x_ii, z_ii] = self.plane_model[x_ii, z_ii] * (1 + new_value/100)
                    else:
                        self.plane_model[x_ii, z_ii] = new_value



    def find_closest_index_in_array(self, matrix, target_value):
        diff = np.abs(matrix - target_value)
        index = np.unravel_index(np.argmin(diff), diff.shape)
        return index[0]

    def set_up_fault_zone(self, surface_x, dip_angle, deepest_z, app_thick, specified_value, percentage_flag):
        """
        Set up a fault zone with the given parameters.
        Args:
            surface_x (float): x coordinate of the intersection between free surface and fault line
            dip_angle (float): dip angle of the fault in degrees
            deepest_z (float): deepest z coordinate of the fault
            app_thick (float): apparent thickness of the fault zone
            specified_value (float): value to set in the fault zone
            percentage_flag (bool): if True, set the specifed value as a percentage of the original velocity
                                    True: NEW = OLD * (1 + specified_value/100)
                                    False: NEW = specified_value
        """
        slope = np.tan(np.deg2rad(dip_angle))

        x_arr = np.arange(self.x_min, self.x_max, self.grid_size) + self.grid_size/2
        z_arr = np.arange(self.z_min, self.z_max, self.grid_size) + self.grid_size/2
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
                    if percentage_flag:
                        self.plane_model[ii, root_min_ii:root_max_ii] = self.plane_model[ii, root_min_ii:root_max_ii] * (1+specified_value/100)
                    else:    
                        self.plane_model[ii, root_min_ii:root_max_ii] = specified_value
    

    def set_up_irregular_zone(self, points, new_value, percentage_flag, s, k, smooth_sigma=3, num_interpolated_points=100):
        """
        Add an irregular zone to the model using spline interpolation to create a boundary with adjustable smoothness.
        
        Args:
            arr (np.array): 2D velocity model array
            points (list of tuples): [(x1, z1), (x2, z2), ..., (xn, zn)] control points defining the region
            new_value (float): New value to set (or percentage change)
            percentage_flag (bool): True for percentage change, False for fixed value
            smooth_sigma (float): Gaussian smoothing parameter to control boundary smoothness
            num_interpolated_points (int): Number of interpolated points for a smoother boundary
            spline_smoothness (float): Spline smoothness (s parameter), 0 for exact fit, larger values for smoother curves
        """
        
        # Convert (x, z) coordinates to index coordinates
        x_list, z_list = zip(*points)
        x_indices = ((np.array(x_list) - self.x_min) // self.grid_size).astype(int)
        z_indices = ((np.array(z_list) - self.z_min) // self.grid_size).astype(int)
        
        # Use B-spline interpolation to smooth the boundary with adjustable smoothness
        tck, u = splprep([x_indices, z_indices], s=s, k=k, per=True)
        u_fine = np.linspace(0, 1, num_interpolated_points)
        x_smooth, z_smooth = splev(u_fine, tck)
        
        # Create region mask
        polygon_path = Path(np.column_stack([x_smooth, z_smooth]))
        mask = np.zeros_like(self.plane_model, dtype=bool)
        
        for x_ii in range(self.nx):
            for z_ii in range(self.nz):
                if polygon_path.contains_point((x_ii, z_ii)):
                    mask[x_ii, z_ii] = True
        
        # Smooth the mask boundary
        mask_smooth = gaussian_filter(mask.astype(float), sigma=smooth_sigma) > 0.5
        
        # Apply new values
        if percentage_flag:
            self.plane_model[mask_smooth] *= (1 + new_value / 100)
        else:
            self.plane_model[mask_smooth] = new_value
    
    def blur_model(self, mode='nearest', sigma=5):
        """
        Apply Gaussian smoothing to the model.
        
        Args:
            mode (str): Mode for Gaussian filter. Default is 'nearest'.
            sigma (float): Standard deviation for Gaussian kernel.
        """
        self.plane_model = gaussian_filter(self.plane_model, mode=mode, sigma=sigma)

    def finalize_model(self):
        """
        Now we need to pad the model with 3 grids each side because of the setting in OpenSWPC
        """
        self.plane_model = np.pad(self.plane_model, ((3, 3), (3, 3)), mode='edge')
        self.nx += 6
        self.nz += 6

    def output_2d(self):
        f = open(f'model_Vp_2d.dat', 'w+')
        print(self.nx, self.nz)
        for ii in range(self.nx):
            for kk in range(self.nz):
                f.write(f'{self.plane_model[ii, kk]:.3f}\n')
        f.close()

    def visualizing_array(self, vmin, vmax, label_name):
        fig, ax = plt.subplots(1,1, figsize=(5,5))
        im = plt.imshow(self.plane_model[:,:].T, cmap='coolwarm_r', extent=[self.x_min, self.x_max, self.z_max, self.z_min], vmin=vmin, vmax=vmax)
        plt.xlabel('X(km)')
        plt.ylabel('Z(km)')
        cbar = plt.colorbar(im, shrink=.4)
        cbar.set_label(label_name)
        fig.savefig('XZ.png', dpi=300)