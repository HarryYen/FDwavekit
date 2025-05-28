# FDwavekit: Simple Note of Creating a 2D Synthetic Model

Here we show you how to create a simple 2D model which can be used when conducting 2D waveform simulation in OpenSWPC.

We can only modify **model_main_2d.py** to create your model, and this script will use the class in synmodel_creator.py to do the work.

In model_main_2d.py: 

1. **setup the model parameters (please keep the same with the setting of your OpenSWPC)**
    
    ```python
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
    ```
    
    Here we specify a 1D velocity model file ‘ak135.csv’, and then your model will be constructed on the basis of this 1D model.
    
    the format of ak135.csv will be:
    
    ```python
    -99 1.45
    0. 1.45
    3. 1.65
    3.3 5.8
    10. 6.8
    18. 8.0355
    43. 8.0379
    80. 8.0450
    120. 8.0505
    165. 8.1750
    ```
    
    the 1st column is depth, and the second one is velocity.
    
    You can also use your own 1D model file but remember to keep the format the same.
    
    Now, your model will be like this:
    
    ![image](https://github.com/HarryYen/FDwavekit/blob/main/images/image.png)
    
2. **put some anomalies in your model**
    
    (1) A normal rectangle
    
    ```python
       synmodel_creator.set_up_rectangle(
            x1 = -20, z1 = 20, x2 = 20, z2 = 20, x3 = 20, z3 = 40, x4 = -20, z4 = 40,  
            new_value = 7., percentage_flag = False
        )
    ```
    
    x1-x4, z1-z2: the position of your 4 points specifying the rectangle.
    
    Result:
    
    ![image.png](https://github.com/HarryYen/FDwavekit/blob/main/images/image%201.png)
    
    Please note that you can specify a perturbation if the percentage_flag is True. 
    
    ```python
        synmodel_creator.set_up_rectangle(
            x1 = -20, z1 = 20, x2 = 20, z2 = 20, x3 = 20, z3 = 40, x4 = -20, z4 = 40,  
            new_value = 7., percentage_flag = True
        )
    ```
    
    Now, the value will be the original value * (1 + 7%).
    
    Result:
    
    ![image.png](https://github.com/HarryYen/FDwavekit/blob/main/images/image%202.png)
    
    (2) Fault zone or oblique layer
    
    ```python
        synmodel_creator.set_up_fault_zone(surface_x = -20,
                                           dip_angle = 20, 
                                           deepest_z = 40., 
                                           app_thick = 10., 
                                           specified_value = 8., 
                                           percentage_flag = False)
    ```
    
    Results:
    
    ![image.png](https://github.com/HarryYen/FDwavekit/blob/main/images/image%203.png)
    

The way we create this oblique layer is like:

![tutorial.png](https://github.com/HarryYen/FDwavekit/blob/main/images/tutorial.png)

Finally you will get a output file called ‘model_Vp_2d.dat’.
