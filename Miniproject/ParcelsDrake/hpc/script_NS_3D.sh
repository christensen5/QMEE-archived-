#!/bin/bash
# Job name
#PBS -N NS_3D_test
# Time required in hh:mm:ss
#PBS -l walltime=00:10:00
# Resource requirements
#PBS -l select=1:ncpus=4:mpiprocs=4:ompthreads=1:mem=15999Mb
# Files to contain standard error and standard output
#PBS -o stdout
#PBS -e stderr
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

# Start time
echo Start time is `date` > data/time

mpiexec python NavierStokes/HPCscript_H_3D.py

# End time
echo End time is `date` >> data/time

# Copy results to WORKDIR
cp -r $TMPDIR/data/* $PBS_O_WORKDIR
