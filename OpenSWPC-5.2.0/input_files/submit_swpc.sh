#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N swpc
#PBS -j oe
#PBS -o mpi.log

###########################################################
# USER PARAMETERS FOR PARALLELIZATION
#PBS -l nodes=node01:ppn=20

##PBS -l nodes=node01:ppn=24+node02:ppn=24+node03:ppn=24+node04:ppn=24+node05:ppn=24+node06:ppn=24+node07:ppn=24+node08:ppn=24+node09:ppn=24+node10:ppn=24+node11:ppn=24
###########################################################


file_list=("test01")

for file in ${file_list[@]};
do
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
exe=../bin/swpc_3d.x
inf=./$file/input.inf
mpirun -np 20 --bind-to core ${exe} -i ${inf}
#mpirun -np 384  ${exe} -i ${inf}
done
