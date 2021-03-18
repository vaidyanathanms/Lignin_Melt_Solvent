#!/bin/bash

#BSUB -P BIP189
#BSUB -W 01:50
#BSUB -nnodes 1
#BSUB -J py_jobname
#BSUB -o outdir/out.%J
#BSUB -e outdir/err.%J



module load caascade/1.1.beta .gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121 gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121
module load gromacs/2020.2-rdtscp_off

export GMX_MAXBACKUP=-1;
export OMP_NUM_THREADS=7;

echo "begin job.."
echo $PWD

mkdir -p outdir
mkdir -p trajfiles

#---------------------------------------------------------Generate index files----------------------------------
ftc_grp=./tcgrp_indx.ndx
if ! test -f "$ftc_grp"; then
	echo "begin generating tempearture coupling groups.."
	# generate temperature_coupling files
	jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi select -s py_finconf -sf py_indexfyle -on tcgrp_indx.ndx
	
fi
#---------------------------------------------------------Generate enermin files----------------------------------
finit_inp=./enermin.tpr
if ! test -f "$finit_inp"; then
	echo "begin generating enermin.tpr.."
	# generate enermin files
	jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi grompp -f minim.mdp -c py_finconf -p py_topol -o enermin.tpr
	wait

fi
wait

#----------------------------------------------------------Minimization--------------------------------------------
fmin_inp=./nvt.tpr 
if ! test -f "$fmin_inp"; then

	echo "begin running enermin.tpr.."
	# run enermin.tpr
	jsrun -X 1 -n 1 -c 42 -a 6 -g 6 --launch_distribution plane:6 -b packed:7 gmx_mpi mdrun -s enermin.tpr -cpo state_min.cpt -cpi state_min.cpt -cpt 2 -g md_min.log -o traj_min.trr -e ener_min.edr -c confout_min.gro -maxh 1.75
	wait

	echo "begin generating nvt.tpr.."
	# generate nvt files
	jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi grompp -f nvt.mdp -c confout_min.gro -p py_topol -n tcgrp_indx.ndx -o nvt.tpr 
	wait

        cp md_min.log trajfiles/md_min.log
	cp traj_min.trr trajfiles/traj_min.trr
	cp ener_min.edr trajfiles/ener_min.edr
	cp confout_min.gro trajfiles/confout_min.gro
fi
wait

#----------------------------------------------------------NVT Equilibration--------------------------------------------
fnvt_inp=./npt_berendsen.tpr
if ! test -f "$fnvt_inp"; then

	echo "begin running nvt.tpr.."
	# run nvt.tpr
	jsrun -X 1 -n 1 -c 42 -a 6 -g 6 --launch_distribution plane:6 -b packed:7 gmx_mpi mdrun -s nvt.tpr -cpo state_nvt.cpt -cpi state_nvt.cpt -cpt 5 -g md_nvt.log -o traj_nvt.trr -e ener_nvt.edr -c confout_nvt.gro -pme gpu -npme 1 -nb gpu -bonded gpu -pin off -maxh 1.75
	wait

	echo "begin generating npt_berendsen.tpr.."
	# generate npt_berendsen files
	jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi grompp -f npt_berendsen.mdp -c confout_nvt.gro -p py_topol -n tcgrp_indx.ndx -o npt_berendsen.tpr
	wait

	cp md_nvt.log trajfiles/md_nvt.log
	cp traj_nvt.trr trajfiles/traj_nvt.trr
	cp ener_nvt.edr trajfiles/ener_nvt.edr
	cp confout_nvt.gro trajfiles/confout_nvt.gro
fi
wait

#----------------------------------------------------------NPT Equilibration--------------------------------------------
fnpt_inp=./npt_main.tpr
if ! test -f "$fnpt_inp"; then

	echo "begin running npt_berendsen.tpr.."
	# run npt_berendsen.tpr
	jsrun -X 1 -n 1 -c 42 -a 6 -g 6 --launch_distribution plane:6 -b packed:7 gmx_mpi mdrun -s npt_berendsen.tpr -cpo state_npt_berendsen.cpt -cpi state_npt_berendsen.cpt -cpt 5 -g md_npt_berendsen.log -o traj_npt_berendsen.trr -e ener_npt_berendsen.edr -c confout_npt_berendsen.gro -pme gpu -npme 1 -nb gpu -bonded gpu -pin off -maxh 1.75
	wait

	echo "begin generating npt_main.tpr.."
	# generate npt_berendsen files
	jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi grompp -f npt_main.mdp -c confout_npt_berendsen.gro -p py_topol -n tcgrp_indx.ndx -o npt_main.tpr
	wait

	cp md_npt_berendsen.log trajfiles/md_npt_berendsen.log
	cp traj_npt_berendsen.trr trajfiles/traj_npt_berendsen.trr
	cp ener_npt_berendsen.edr trajfiles/ener_npt_berendsen.edr
	cp confout_npt_berendsen.gro trajfiles/confout_npt_berendsen.gro
else

        echo "begin running npt_main.tpr.."
        # run npt_main.tpr
        jsrun -X 1 -n 1 -c 42 -a 6 -g 6 --launch_distribution plane:6 -b packed:7 gmx_mpi mdrun -s npt_main.tpr -cpo state_npt_main.cpt -cpi state_npt_main.cpt -cpt 5 -g md_npt_main.log -o traj_npt_main.trr -e ener_npt_main.edr -c confout_npt_main.gro -pme gpu -npme 1 -nb gpu -bonded gpu -pin off -maxh 1.75
        wait


        cp md_npt_main.log trajfiles/md_npt_main.log
        cp traj_npt_main.trr trajfiles/traj_npt_main.trr
        cp ener_npt_main.edr trajfiles/ener_npt_main.edr
        cp confout_npt_main.gro trajfiles/confout_npt_main.gro
fi
wait
#-----------------------------------High Temperature NVT----------------------------------------------------------------------------

fnvthigh_inp=./nvt_high.tpr 
if ! test -f "$fnvthigh_inp"; then

	echo "begin generating high temperature NVT: nvt_high.tpr.."
	# generate nvt_high files
	jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi grompp -f nvt_high.mdp -c confout_npt_main.gro -p py_topol -n tcgrp_indx.ndx -o nvt_high.tpr
	wait

	echo "begin running high temperature NVT: nvt_high.tpr.."
	# run nvt_high.tpr
	jsrun -X 1 -n 1 -c 42 -a 6 -g 6 --launch_distribution plane:6 -b packed:7 gmx_mpi mdrun -s nvt_high.tpr -cpo state_nvt_high.cpt -cpi state_nvt_high.cpt -cpt 5 -g md_nvt_high.log -o traj_nvt_high.trr -e ener_nvt_high.edr -c confout_nvt_high.gro -pme gpu -npme 1 -nb gpu -bonded gpu -pin off -maxh 1.75
	wait


	cp md_nvt_high.log trajfiles/md_nvt_high.log
	cp traj_nvt_high.trr trajfiles/traj_nvt_high.trr
	cp ener_nvt_high.edr trajfiles/ener_nvt_high.edr
	cp confout_nvt_high.gro trajfiles/confout_nvt_high.gro
else

        echo "begin running high temperature NVT: nvt_high.tpr.."
        # run nvt_high.tpr
        jsrun -X 1 -n 1 -c 42 -a 6 -g 6 --launch_distribution plane:6 -b packed:7 gmx_mpi mdrun -s nvt_high.tpr -cpo state_nvt_high.cpt -cpi state_nvt_high.cpt -cpt 5 -g md_nvt_high.log -o traj_nvt_high.trr -e ener_nvt_high.edr -c confout_nvt_high.gro -pme gpu -npme 1 -nb gpu -bonded gpu -pin off -maxh 1.75
        wait


        cp md_nvt_high.log trajfiles/md_nvt_high.log
        cp traj_nvt_high.trr trajfiles/traj_nvt_high.trr
        cp ener_nvt_high.edr trajfiles/ener_nvt_high.edr
        cp confout_nvt_high.gro trajfiles/confout_nvt_high.gro
fi
wait
#------------------------------------------------------------------------------------------------------------------------------------

echo "End of run.."

