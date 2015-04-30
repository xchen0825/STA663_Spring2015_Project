#!/bin/bash
# parallel job using -n processors. and runs for 72 hours (max)
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=72000
#SBATCH --ntasks=2
#SBATCH --job-name=mpi_test
#SBATCH --output=mpi_test.out
#SBATCH --error=mpi_test.out
export I_MPI_PMI_LIBRARY=/opt/slurm/lib64/libpmi.so
source /opt/apps/intel/intelvars.sh
srun -n $SLURM_NTASKS inspxe-runmc /dscrhome/xc46/PIHM_IntelMPI/pihm -n 2 
