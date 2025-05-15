from model_main import nx, ny, nz, grid_size
from model_main import x_min, z_min, x_max, z_max, y_min, y_max
import numpy as np
from matplotlib.path import Path
import matplotlib.pyplot as plt

def find_closest_index_in_array(matrix, target_value):
    diff = np.abs(matrix - target_value)
    index = np.unravel_index(np.argmin(diff), diff.shape)
    return index[0]

def set_up_fault_plane(fault_x, fault_dip, fault_Vp, model_arr):
    """
    Usage:
        It is used to create a fault plane in a 2D model.
    Args:
        fault_x   : the location of the fault on the surface.
        fault_dip : the dip of the fault.
        fault_Vp  : the velocity of P wave of this fault.
        model_arr : the initial 2-D model array
    Returns:
        model_arr : the updated model array.
    """
    global nx, nz, x_min, z_min, x_max, z_max, grid_size
    slope = np.tan(np.deg2rad(fault_dip))

    def fault_equation(x):
        return slope * (-x + fault_x)

    x_arr = np.arange(x_min, x_max, grid_size) + grid_size/2
    z_arr = np.arange(z_min, z_max, grid_size) + grid_size/2
    root_arr = -1 * fault_equation(x_arr)

    for x_ii, root_z in enumerate(root_arr):
        if root_z < z_min or root_z > z_max:
            continue
        z_ii = find_closest_index_in_array(z_arr, root_z)
        model_arr[x_ii, z_ii] = fault_Vp
    
    return model_arr


def set_up_fault_zone(fault_x, fault_dip, fault_Vp, app_thick, model_arr, deepest_z):
    """
    Usage:
        It is used to create a fault plane in a 2D model.
    Args:
        fault_x   : the location of the fault on the surface.
        fault_dip : the dip of the fault.
        fault_Vp  : the velocity of P wave of this fault.
        model_arr : the initial 2-D model array
        deepest_z  : the deepest depth of the fault (km)
    Returns:
        model_arr : the updated model array.
    """
    global nx, nz, x_min, z_min, x_max, z_max, grid_size
    slope = np.tan(np.deg2rad(fault_dip))

    def fault_eq1(x):
        return slope * (-x + fault_x)
    
    def fault_eq2(x):
        return slope * (-x + fault_x) + app_thick

    x_arr = np.arange(x_min, x_max, grid_size) + grid_size/2
    z_arr = np.arange(z_min, z_max, grid_size) + grid_size/2
    z_arr_index = range(nz)
    root1_arr = -1 * fault_eq1(x_arr)
    root2_arr = -1 * fault_eq2(x_arr)

    for ii in range(len(x_arr)):
        root_min = min(root1_arr[ii], root2_arr[ii])
        root_max = max(root1_arr[ii], root2_arr[ii])
        if root_max > deepest_z:
            if root_min > deepest_z: 
                continue
            elif root_min <= deepest_z:
                root_max = deepest_z
        
        if not ((root_min <= z_min and root_max <= z_min) or (root_min >= z_max and root_max >= z_max)):
            root_min_ii = find_closest_index_in_array(z_arr, root_min)
            root_max_ii = find_closest_index_in_array(z_arr, root_max)
            
            if (root_min_ii == z_arr_index[0] and root_max_ii == z_arr_index[0]) or (root_min_ii == z_arr_index[-1] and root_max_ii == z_arr_index[-1]):
                continue
            else:    
                model_arr[ii, root_min_ii:root_max_ii] = fault_Vp
    return model_arr


# def is_inside_rectangle(x, y, x1, y1, x2, y2, x3, y3, x4, y4):
#     # Check if point (x, y) is inside the rectangle defined by four linear lines.
#     def is_left_of_line(x, y, x1, y1, x2, y2):
#         return (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) >= 0

#     return is_left_of_line(x, y, x1, y1, x2, y2) and is_left_of_line(x, y, x2, y2, x3, y3) \
#         and is_left_of_line(x, y, x3, y3, x4, y4) and is_left_of_line(x, y, x4, y4, x1, y1)

# def set_up_a_rectangle(arr, x1, y1, x2, y2, x3, y3, x4, y4, grid_size, new_value):
#     # Iterate over the 2D array and update values inside the specified rectangle.
#     x_list = np.array([x1, x2, x3, x4])
#     y_list = np.array([y1, y2, y3, y4])
#     x_ii_list = (x_list - x_min) // grid_size + 1
#     y_ii_list = (y_list - z_min) // grid_size + 1
#     x1_ii, x2_ii, x3_ii, x4_ii = x_ii_list[0], x_ii_list[1], x_ii_list[2], x_ii_list[3]
#     y1_ii, y2_ii, y3_ii, y4_ii = y_ii_list[0], y_ii_list[1], y_ii_list[2], y_ii_list[3]
#     print(x_ii_list)
#     print(y_ii_list)
#     for x in range(nx):
#         for y in range(nz):
#             if is_inside_rectangle(x, y, x1_ii, y1_ii, x2_ii, y2_ii, x3_ii, y3_ii, x4_ii, y4_ii):
#                 print(x,y)
#                 arr[x,y] = new_value
#     return arr

def set_up_rectangle(arr, x1, z1, x2, z2, x3, z3, x4, z4, new_value):
    x_list = np.array([x1, x2, x3, x4])
    z_list = np.array([z1, z2, z3, z4])
    x_ii_list = (x_list - x_min) // grid_size + 1
    z_ii_list = (z_list - z_min) // grid_size + 1
    x1_ii, x2_ii, x3_ii, x4_ii = x_ii_list[0], x_ii_list[1], x_ii_list[2], x_ii_list[3]
    z1_ii, z2_ii, z3_ii, z4_ii = z_ii_list[0], z_ii_list[1], z_ii_list[2], z_ii_list[3]
    p = Path([(x1_ii, z1_ii), (x2_ii, z2_ii), (x3_ii, z3_ii), (x4_ii, z4_ii)])
    for x_ii in range(nx):
        for z_ii in range(nz):
            is_in_rectangle = p.contains_point((x_ii, z_ii))
            if is_in_rectangle:
                print(x_ii, z_ii)
                arr[x_ii, z_ii] = new_value
    return arr

def plot_2D_figures(arr, x, y, z):
    
    x_ii = int((x - x_min) // grid_size) + 4  # extra 3 points for padding
    y_ii = int((y - y_min) // grid_size) + 4
    z_ii = int((z - z_min) // grid_size) + 4
    comp_list = ['x', 'y', 'z']

    fig, ax = plt.subplots(1,1, figsize=(5,5))
    plt.imshow(arr[x_ii, :, :].T, cmap='RdYlBu', extent=(y_min,y_max,z_max,z_min), vmin=1.5, vmax=3)
    plt.xlabel('Y(km)')
    plt.ylabel('Z(km)')
    plt.colorbar(shrink=.4)
    plt.savefig('YZ.png', dpi=300)
    plt.close()
    
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    plt.imshow(arr[:, y_ii, :].T, cmap='RdYlBu', extent=(x_min,x_max,z_max,z_min), vmin=1.5, vmax=3)
    plt.xlabel('X(km)')
    plt.ylabel('Z(km)')
    plt.colorbar(shrink=.4)
    plt.savefig('XZ.png', dpi=300)
    plt.close()
    
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    plt.imshow(arr[:, :, z_ii].T, cmap='RdYlBu', aspect=5, extent=(x_min,x_max,y_min,y_max), vmin=1.5, vmax=3)
    plt.xlabel('X(km)')
    plt.ylabel('Y(km)')
    plt.colorbar(shrink=.4)
    plt.savefig('XY.png', dpi=300)
    plt.close()
    
    
    
    


def output(Vel_arr):
    f = open(f'model_Vp.dat', 'w+')
    for jj in range(ny + 6):
        for ii in range(nx + 6):
            for kk in range(nz + 6):
                f.write(f'{Vel_arr[ii, jj, kk]:.3f}\n')
    f.close()
    
def output_2d(Vel_arr):
    f = open(f'model_Vp_2d.dat', 'w+')
    print(nx, nz)
    for ii in range(nx + 6):
        for kk in range(nz + 6):
            f.write(f'{Vel_arr[ii, kk]:.3f}\n')
    f.close()