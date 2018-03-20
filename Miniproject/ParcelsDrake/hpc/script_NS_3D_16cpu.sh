#!/bin/bash
# Job name
#PBS -N NS_3D_16cpu
# Time required in hh:mm:ss
#PBS -l walltime=24:00:00
# Resource requirements
#PBS -l select=1:ncpus=16:mem=15999Mb
# Files to contain standard error and standard output
#PBS -o stdout16
#PBS -e stderr16
# Mail notification
#PBS -m ae
#PBS -M akc17@imperial.ac.uk

echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
rm -f stdout* stderr*

module load gcc
module load mpi
export I_MPI_CC=gcc
export I_MPI_CXX=g++

# Change the next line to wherever you put your firedrake/
source $HOME/firedrake/bin/activate

# Switch to TMPDIR
cp -r $PBS_O_WORKDIR/* $TMPDIR
cd $TMPDIR
mkdir -p data

# Start time
echo Start time is `date` > data/time16

mpiexec -n 16 python NavierStokes/HPCscript_H_3D.py

# End time
echo End time is `date` >> data/time16

# Copy results to WORKDIR
cp -r data $PBS_O_WORKDIR
