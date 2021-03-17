#!/bin/bash

#BSUB -P BIP189
#BSUB -W 00:10
#BSUB -nnodes 1
#BSUB -J py_jobname
#BSUB -o outdir/out.%J
#BSUB -e outdir/err.%J

module load caascade/1.1.beta .gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121 gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121
module load gromacs/2020.2-rdtscp_off

export GMX_MAXBACKUP=-1;
export OMP_NUM_THREADS=2;

echo "begin job.."
echo $PWD

mkdir -p initdir

# editconf box
jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi editconf -f py_meltconf -bt cubic -d 0.5  -o py_boxmeltconf
wait

# solvate with solvent - organic solvent
py_solvate_1

# solvate with cosolvent - water
py_solvate_2

# make enermin_tpr file
jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi grompp -f minim.mdp -p py_topol -c py_finconf -o enermin.tpr
wait

# generate temperature_coupling files
jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi select -s py_finconf -sf py_indexfyle -on tcgrp_indx.ndx
wait

cp *.pdb initdir/
cp *.psf initdir/
cp *.top initdir/
cp *.txt initdir/
cp *.ndx initdir/
cp *.gro initdir/
