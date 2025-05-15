#!/bin/bash

cd swpc_3d/
make clean
cp ../netcdf.mod .
#cp netcdf.mod swpc_sh
#cp netcdf.mod swpc_3d
#cp netcdf.mod tools
make arch=ubuntu-gfortran NCLIB=-L/cluster/gcc630/netcdf-4.4.1.1-nc4/lib
#cp ../interp_test/linear_interp.mod .

#cd ../../
#mpirun -np 2 ./bin/swpc_psv.x -i input.inf
#cd src
