# Interpolating Models for FDM Waveform Modeling
This md file is written by **Hung-Yu Yen.**

It's a simple tutorial for script MODEL_PREPARATION.

This script will read the velocity (or density) model you prepare and interpolate it into the FDM grids you need in Waveform Modeling.

Here, this file will guide you step-by-step in the following part.

###### Last updating: 2024/09/02, Hung-Yu Yen 
## Some Must-Know

### Why do we need this script?
Tomography models usually have their own grid spacing (e.g. The grid spacing of the P/S Joint Inversion of entire Taiwan Region (Hsin-Hua Huang et al., 2014) is roughly 0.8 km), which might be different from our spacing in waveform modeling.

Therefore, we need interpolation.

## 0. Install the module
```
pip install -r requirements.txt
```
If there are some errors, you can manually check the modules in requirements.txt and install them.

## 0. Preparing your own velocity model
In this script, you need to prepare the velocity file like this:

Please see the example below:
###### model_input/TW-PS-H14.txt
```
Dimension: 76 61 27

lon lat dep vp ....

119.000 21.300 -5.000 3.000 ...
```
***The first row records the dimension of the whole model (the number of grids on longitude axis is 76, latitute for 61 and depth for 27)***

It means that the total grid number of this velocity model is 76 x 61 x 27 = 125,172

***The second row is the header***

**lon lat dep** are required!!!

It doesn't matter that this table contains other header because our script will not use it. **However, you need to confirm the vector you want to interpolate (like vp, vs or rho) is also required.**


## 1. Modify the inpur file
###### input.yaml
```
# input.yaml
dx: 0.15
dy: 0.15
dz: 0.15
cx: 121.7
cy: 24.735
x_half_len: 30.0
y_half_len: 30.0
z_beg: -2.1
z_end: 100.0
interpolation_method: 'linear'
vector: 'vs'
default_value: 1.5
cores_num: 10
```
**dx, dy, dz**: grid spacing in FDM

**cx, cy**: the center (lon, lat) in your computational region in FDM.

**x_half_len, y_half_len**: the **half length** of your model.
```
+-------------------+
|                   |
|                   |
|        *--x_half--|
|                   |
|                   |
+-------------------+


*: the position of (cx, cy)
```
**z_beg**: the shallowest depth in your computational region. (**It is required to be negative in OpenSWPC!!!**)

**z_end**: the deepest depth in your computational region.

**interpolation_method**: linear, nearest

**vector**: the vector you want to interpolate (e.g. 'vp', 'vs', 'rho'). (**It need to be the same as the header name in your input model**)

**default_value**: this value will be filled if the computational region you set is out of the input model.

**cores_num**: cores number for multiprpcessing.

## 2. executing
```
python main.py
```

## 3. output file
Output file will be in the same directory of main.py.
```
./3d_model_{vector_name}.dat
```

### Something you need to know
1. the interpolated file (3d_model_{vector_name}.dat) has slightly bigger range than the range you specify (x_half_len, y_half_len, cx, cy), because OpenSWPC script adds 3 grids more in each direction.

    Says, if dx is 0.1 and x_half_len is 5, we should have 5 x 2 / 0.1 = 100 grids on x-direction. **But** OpenSWPC will add 3 grids in both before start point and after end point. Finally, we got 100 + 3 + 3 = 106 grids. 
    
    Therefore, the model range has 53 grids before cx and also 53 grids after cx.

