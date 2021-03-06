# To analyze lignin results
# Lignin-melt and lignin-solvent systems. 
# Version: Mar-20-2021
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
from make_ana import * # function definitions
#------------------------------------------------------------------

# Version Info and # of input args for parsing thermostat coeff
print("Generating GROMACS analysis codes")
print("Version: Mar-20-2021")
#------------------------------------------------------------------

# Input Data
run_ana   = 1 # 1-copy files and run, 0-NO run (copies files)
inp_type  = 'cosolvents' # melts, solvents, cosolvents
biomass   = 'MYB' # name of the biomass type
disperse  = 'mono' # mono/poly; only for melts
run_arr   = [4,5,6] # number of independent runs for a given biomass
otyp_arr  = ['EOH','THF','GVL']  # solvent arr for solvents/cosolvents
oname_arr = otyp_arr # change if prefix is different from name in PDB
wtyp_arr  = ['tip3p','tip3p','tip3p'] # water arr type
wname_arr = ['TIP3_','TIP3_','TIP3_'] # diff from prefix
nchains   = 1     # number of polymer chains
npoly_res = 22  # number of polymer residues
#------------------------------------------------------------------

# Directory Paths
main_dir  = os.getcwd() # current dir
gmx_dir   = '../src_gmx' # gmx file super directory
ana_dir   = gmx_dir + '/' + 'ana_files' # sh file dir
scr_dir   = '/gpfs/alpine/bip189/scratch/vaidyams' # scratch dir

if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
scr_dir  = scr_dir + '/lignin'
if not os.path.isdir(scr_dir):
    os.mkdir(scr_dir)
#------------------------------------------------------------------

# Required GMX/sh Files
sh_anafyl = 'runana_pyinp.sh'
ana_fyles = ['nneigh_inp.txt','rdfinp1.txt','rdfinp2.txt',\
             'rg_inp.txt','hbinp.txt','sasainp.txt']
#------------------------------------------------------------------

# Check input dimension consistency
if inp_type == 'cosolvents':
    check_arr_dim(len(otyp_arr),len(oname_arr),len(wtyp_arr),\
                  len(wname_arr))
    solv_len = len(otyp_arr)
elif inp_type == 'melts':
    solv_len = 0
#------------------------------------------------------------------

#Main Code
for inp_val in range(solv_len): # loop in solvents

    o_sol_typ,solv_name,wat_type,wat_name=assign_vals(otyp_arr,oname_arr,\
                                                      wtyp_arr,wname_arr,\
                                                      inp_type,inp_val)

    for casenum in range(len(run_arr)): # loop in runarr

        print( "Analyzing: ", biomass,inp_type,o_sol_typ,\
               run_arr[casenum])

        # Make directories
        head_dir = scr_dir + '/' + inp_type
        if not os.path.isdir(head_dir):
            print(head_dir, " does not exist")
            continue

        poly_dir = head_dir + '/' + biomass
        if not os.path.isdir(poly_dir):
            print(poly_dir, " does not exist")
            continue

        if inp_type == 'melts':
            poly_dir = poly_dir + '/' + disperse
            if not os.path.isdir(poly_dir):
                print(poly_dir, " does not exist")
                continue

        rundir = poly_dir + '/run_' + str(run_arr[casenum])
        if not os.path.isdir(rundir):
            print(rundir, " does not exist")
            continue

        workdir1 = set_working_dir(rundir,inp_type,o_sol_typ)

        # Check and copy ana_files/sh_file
        check_cpy_ana_files(ana_fyles,sh_anafyl,inp_type,o_sol_typ\
                            ,ana_dir,workdir1)
        
        # retrieve tpr/trr file NEED TO WORK ON THIS
        #        tpr_file,trr_file=retrieve_traj_files(workdir1)


        os.chdir(workdir1)

        # Make an outdir for writing output files from summit
        sum_out_dir = workdir1 + '/outdir'
        if not os.path.isdir(sum_out_dir):
            os.mkdir(sum_ana_dir)

        # Edit runana shell script file
        edit_anash_file(biomass,inp_type,o_sol_typ,sh_anafyl)

        # Submit preprocess job
        if run_ana:
            print("Submitting analysis: runana.sh")
            subprocess.call(["chmod", "777", "runana.sh"])
            subprocess.call(["bsub", "runana.sh"])

        # Write end of loop and return directory handle to main dir
        print( "Completed run ", run_arr[casenum])
        os.chdir(main_dir)# current dir
