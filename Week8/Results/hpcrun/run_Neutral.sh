#!/bin/bash
#PBS -lselect=1:ncpus=1:mem=1gb
#PBS -lwalltime=12:00:00
module load R
module load intel-suite
echo "R is about to run"
R --vanilla < $WORK/model_functions_HPC.R
mv AKC_NTSresult_iter* $WORK
echo "R has finished running"
