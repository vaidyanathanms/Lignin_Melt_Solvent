# To run different initial conditions for polydisperse/monodisperse
# lignin and lignin-solvent systems. 
# Version: Mar-03-2021
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
from make_gmx_inputs import * # function definitions
#------------------------------------------------------------------

# Version Info and number of input arguments
print("Generating GROMACS run-time inputs")
print("Version: Mar-03-2021")
print 
#------------------------------------------------------------------

# Input Data
inp_type  = 'melts' # melts, solvents, cosolvents
biomass   = 'switchgrass' # name of the biomass type
disperse  = 'mono'  # mono/poly; only for melts
num_cases = 5 # for solvents, this will be == len(norg_solv)
nchains   = 50 # number of chains in the polymer 
norg_solv = [1000] # number of organic solvents
nwater    = [4000] # number of water molecules (for cosolvents)
#------------------------------------------------------------------

# Directory Paths
main_dir = os.getcwd()
gmx_dir  = '../src_gmx'
mdp_dir  = gmx_dir + '/' + 'mdp_files'
sh_dir   = gmx_dir + '/' + 'sh_files'
scr_dir  = '/gpfs/alpine/bip189/scratch/vaidyams/'
if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
scr_dir  = scr_dir + '/' + inp_type
work_dir = '/gpfs/alpine/bip189/scratch/vaidyams/lignin/melts/large_systems'
#------------------------------------------------------------------

# Required GMX Files
mdp_fyles = ['minim_pyinp.mdp','nvt_pyino.mdp',\
             'npt_berendsen_pyinp.mdp','npt_main_pyinp.mdp']
#------------------------------------------------------------------

#Main Code
for casenum in range(num_cases):
    
    print( "Case number: ", casenum + 1)
      
    # Make directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        os.mkdir(head_dir)

    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        os.mkdir(poly_dir)

    disper_dir = poly_dir + '/' + disperse
    if not os.path.isdir(dispdir):
        os.mkdir(disper_dir)
    
    workdir1 = disper_dir + '/casenum_' + str(casenum+1)
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    # Set thermostat coefficients
    temp_vel_resc,temp_par_rah,pres_par_rah = couple_coeff(coeff_fyle,
    \)

    # Change to working directory
    os.chdir(workdir1)

    # Check for pdb/psf/top files
    check_inp_files(workdir1)

    # Add mdp/shell script files
    check_cpy_mdp_files(mdp_dir,workdir1,sh_fyles)
    rerun = cpy_edit_sh_files(sh_dir,workdir1,mdp_fyles)


        
    for iarch in range(len(archarr)):
             
        if archarr[iarch] == 1:
            print( "Archval: Block_Block")
            dirstr = 'bl_bl'
            fylstr = 'block_block'
        elif archarr[iarch] == 2:
            print( "Archval: Block_Alter")
            dirstr = 'bl_al'
            fylstr = 'block_alter'
        elif archarr[iarch] == 3:
            print( "Archval: Alter_Block")
            dirstr = 'al_bl'
            fylstr = 'alter_block'
        elif archarr[iarch] == 4:
            print( "Archval: Alter_Alter")
            dirstr = 'al_al'
            fylstr = 'alter_alter'
        else:
            print( "Unknown Architecture")
            
        workdir_arch = workdir2 + '/' + dirstr

        if not os.path.isdir(workdir_arch):
            print("ERROR: ", workdir_arch, " not found")
            continue

        workdir_freepdi = workdir_arch + '/pdi_free_' + str(pdi_free)
        if not os.path.isdir(workdir_freepdi):
            print("ERROR: ", workdir_freepdi, " not found")
            continue

        workdir_graftpdi = workdir_freepdi + '/pdi_graft_' + str(pdi_graft)
        if not os.path.isdir(workdir_graftpdi):
            print("ERROR: ", workdir_graftpdi, " not found")
            continue

        for casenum in range(len(ncases_pdi)):

            workdir_subpdi = workdir_graftpdi + '/Case_' + str(ncases_pdi[casenum])
            if not os.path.isdir(workdir_subpdi):
                print("ERROR: ", workdir_subpdi, " not found")
                continue

            os.chdir(workdir_subpdi)
            destdir = os.getcwd()

            print( "Starting analysis for case ", ncases_pdi[casenum], "in ",\
                       free_chains[ifree],dirstr)

            #---Copying files------
            print( "Current Dir ", destdir)
            print( "Copying Files")
                
            for fyllist in range(len(ana_files)):
                cpy_main_files(maindir,destdir,ana_files[fyllist])

            for fyllist in range(len(job_files)):
                cpy_main_files(jobdir,destdir,job_files[fyllist])

            dataname =  find_datafyle(free_chains[ifree],fylstr,\
                                      ncases_pdi[casenum],destdir)
            if dataname == 'ERROR':
                print("ERROR: No restart files found")
                continue

            
            #latest_traj = find_latest_trajfyle(traj_pref,destdir)
            #if latest_traj == 'ERROR':
            #    print("ERROR: No trajectory files found")
            #    continue
            
            compile_anafiles()            
            ntotch = free_chains[ifree] + graft_chains
            traj_arr = glob.glob(traj_pref)
            if traj_arr == []:
                print("ERROR: No trajectory files found")
                continue
            
            for fyllist in range(len(traj_arr)):
                print("Analyzing ", traj_arr[fyllist])
                edit_generate_anainp_files(dataname,traj_arr[fyllist],ntotch,\
                                           free_chains[ifree],graft_chains,\
                                           cutoff_dist,fyllist+1)
                outana = 'jobana_' + str(fyllist+1) + '.sh'
                run_analysis(free_chains[ifree],pdi_free,ncases_pdi[casenum],\
                             fylstr,'jobana_var.sh',outana,fyllist+1)

