#!/bin/bash

#BSUB -P BIP189
#BSUB -W 01:55
#BSUB -nnodes 1
#BSUB -J py_jobname
#BSUB -o outdir/anaout.%J
#BSUB -e outdir/anaerr.%J

module load caascade/1.1.beta .gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121 gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121
module load gromacs/2020.2-rdtscp_off

export GMX_MAXBACKUP=-1;
export OMP_NUM_THREADS=7;

echo "begin job.."
echo $PWD

mkdir -p outdir
mkdir -p anafiles

# make index files
jsrun -X 1 -n 1 -c 1 -a 1 -g 1 --launch_distribution plane:1 -b packed:1 gmx_mpi select -s nvt_high.tpr -on nneigh_inp.ndx -sf nneigh_inp.txt
wait

jsrun -X 1 -n 1 -c 1 -a 1 -g 1 --launch_distribution plane:1 -b packed:1 gmx_mpi select -s nvt_high.tpr -on rg_inp.ndx -sf rg_inp.txt
wait

jsrun -X 1 -n 1 -c 1 -a 1 -g 1 --launch_distribution plane:1 -b packed:1 gmx_mpi select -s nvt_high.tpr -on hbinp.ndx -sf hbinp.txt
wait

# compute RDF
jsrun -X 1 -n 1 -c 14 -a 2 -g 2 --launch_distribution plane:2 -b packed:7 gmx_mpi rdf -f traj_nvt_high.trr -s nvt_high.tpr -sf rdfinp1.txt -b 250000 -e 350000 -o rdfout1.xvg &

jsrun -X 1 -n 1 -c 14 -a 2 -g 2 --launch_distribution plane:2 -b packed:7 gmx_mpi rdf -f traj_nvt_high.trr -s nvt_high.tpr -sf rdfinp2.txt -b 250000 -e 350000 -o rdfout2.xvg &

# compute SASA
jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi sasa -f traj_nvt_high.trr -s nvt_high.tpr -sf sasainp.txt -b 50000 -o sasa.xvg &

# compute radius of gyration
jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi gyrate -f traj_nvt_high.trr -s nvt_high.tpr -n rg_inp.ndx -o rg_chain.xvg -b 50000 &

# compute water-beta_O_4 contacts/solvent-beta_O_4 contacts
echo "0 1 2" | jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi mindist -f traj_nvt_high.trr -s nvt_high.tpr -n nneigh_inp.ndx -on bO4_solv_water.xvg -o genmindist.out -b 5000 -ng 2 &


# compute water-polymer hb
echo "0 2" | jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi hbond -f traj_nvt_high.trr -s nvt_high.tpr -n hbinp.ndx -num wat_hbnum1.xvg -dist wat_hbdist1.xvg -b 50000 -e 150000 &

echo "0 2" | jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi hbond -f traj_nvt_high.trr -s nvt_high.tpr -n hbinp.ndx -num wat_hbnum2.xvg -dist wat_hbdist2.xvg -b 150010 -e 250000 &

echo "0 2" | jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi hbond -f traj_nvt_high.trr -s nvt_high.tpr -n hbinp.ndx -num wat_hbnum3.xvg -dist wat_hbdist3.xvg -b 250010 -e 350000 &

# compute org_solv-polymer hb
echo "0 1" | jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi hbond -f traj_nvt_high.trr -s nvt_high.tpr -n hbinp.ndx -num sol_hbnum1.xvg -dist sol_hbdist1.xvg -b 50000 -e 150000 &

echo "0 1" | jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi hbond -f traj_nvt_high.trr -s nvt_high.tpr -n hbinp.ndx -num sol_hbnum2.xvg -dist sol_hbdist2.xvg -b 150010 -e 250000 &

echo "0 1" | jsrun -X 1 -n 1 -c 7 -a 1 -g 1 --launch_distribution plane:1 -b packed:7 gmx_mpi hbond -f traj_nvt_high.trr -s nvt_high.tpr -n hbinp.ndx -num sol_hbnum3.xvg -dist sol_hbdist3.xvg -b 250010 -e 350000 &

wait

echo "completed all analysis"

cp sol_hb* anafiles/
cp wat_hb* anafiles/
cp hbinp* anafiles/
cp *.ndx anafiles/
cp nneigh* anafiles/
cp rdfinp* anafiles/
cp rdfout* anafiles/
cp sasa* anafiles/
cp bO4* anafiles/
cp genminout* anafiles/
cp rg* anafiles/
