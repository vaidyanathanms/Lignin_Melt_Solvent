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

# Version Info and # of input args for parsing thermostat coeff
print("Generating GROMACS run-time inputs")
print("Version: Mar-03-2021")
if len(sys.argv) == 1:
    coeff_fyle = 'None' 
elif len(sys.argv) == 2:
    coeff_fyle = str(sys.argv[1])
else:
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()
#------------------------------------------------------------------

# Input Data
inp_type  = 'cosolvents'       # melts, solvents, cosolvents
biomass   = 'WT' # name of the biomass type
disperse  = 'mono' # mono/poly; only for melts
o_sol_typ = 'EOH'  # prefix for solvent file for solvents/cosolvents
wat_type  = 'tip3p' # prefix for water file for coslvents
solv_name = o_sol_typ # change if prefix is different from name in PDB
wat_name  = 'TIP3_' # diff from prefix
run_arr   = [1,2,3] # run number
nchains   = 1     # number of polymer chains
npoly_res = 22  # number of polymer residues
n_orgsolv = 2492 # number of organic solvents
nwater    = 4063 # number of water molecules (for cosolvents)
box_dim   = 15 # box size for solvent only. cosolvent=+3
#------------------------------------------------------------------

# Directory Paths
main_dir  = os.getcwd() # current dir
gmx_dir   = '../src_gmx' # gmx file super directory
top_dir   = gmx_dir + '/solv_files/topol' # topology dir
cfg_dir   = gmx_dir + '/solv_files/initguess' # configuration dir
prm_dir   = gmx_dir + '/solv_files/prm_itp' # prm file dir
mdp_dir   = gmx_dir + '/' + 'mdp_files' # mdp file dir
sh_dir    = gmx_dir + '/' + 'sh_files' # sh file dir
scr_dir   = '/gpfs/alpine/bip189/scratch/vaidyams' # scratch dir

if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
scr_dir  = scr_dir + '/lignin'
if not os.path.isdir(scr_dir):
    os.mkdir(scr_dir)
#------------------------------------------------------------------

# Required GMX/sh Files
mdp_fyles  = ['minim_pyinp.mdp','nvt_pyinp.mdp',\
              'npt_berendsen_pyinp.mdp','npt_main_pyinp.mdp']
sh_md_fyle = 'run_md_pyinp.sh'
sh_pp_fyle = 'run_preprocess_pyinp.sh'
sh_rep_fyl = 'repeat.sh'
#------------------------------------------------------------------

#Main Code
for casenum in range(len(run_arr)):
    
    print( "Run number: ", run_arr[casenum])
      
    # Make directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        print(head_dir, " does not exist")
        print("ERR: Create path")
        continue

    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        print(poly_dir, " does not exist")
        print("ERR: Create path")
        continue

    if inp_type == 'melts':
        poly_dir = poly_dir + '/' + disperse
        if not os.path.isdir(poly_dir):
            print(poly_dir, " does not exist")
            print("ERR: Create path")
            continue
    
    rundir = poly_dir + '/run_' + str(run_arr[casenum])
    if not os.path.isdir(rundir):
        print(rundir, " does not exist")
        print('ERR: Create path and input top/pdb files first')
        continue

    workdir1 = set_working_dir(rundir,inp_type,o_sol_typ)

    # Set thermostat/top variables (change if there is a temp sweep)
    Tetau_nvt,Tetau_berend,Tetau_parrah,Prtau_berend,Prtau_parrah,\
        ref_temp,ref_pres,melt_topname = couple_coeff(inp_type,\
                                                      coeff_fyle)

    # Check for pdb/psf/top files of the melt/polymer
    poly_conffile,poly_topfile=check_inp_files(workdir1,melt_topname)

    # Find tc_groups
    tc_grp,tc_typ=create_tcgrps(inp_type,npoly_res,solv_name,wat_name)
    
    # Copy and edit mdp files
    check_cpy_mdp_files(mdp_dir,workdir1,mdp_fyles,inp_type,Tetau_nvt\
                        ,Tetau_berend,Tetau_parrah,Prtau_berend,\
                        Prtau_parrah,ref_temp,ref_pres,tc_grp,tc_typ,\
                        main_dir,coeff_fyle)
    cont_run,edit_sh_fyle = cpy_sh_files(sh_dir,workdir1,\
                                         sh_pp_fyle,sh_md_fyle)


    # Copy/edit top/conf/prm files for solvents/cosolvents
    ff_dir = create_ff_dir(workdir1)
    if inp_type == 'solvents' or inp_type == 'cosolvents':
        if cont_run == 0:
            # sol_top/sol_prm/sol_cfg are arrays
            sol_top,sol_prm,sol_cfg=cpy_solv_files(top_dir,cfg_dir,\
                                                   prm_dir,ff_dir,\
                                                   inp_type,o_sol_typ,\
                                                   wat_type)
            
            poly_topedit = edit_main_top_file(poly_topfile,ff_dir,\
                                              sol_top,sol_prm,workdir1)

    # Copy and edit repeat.sh first time
    if not os.path.exists(workdir1 + '/' + sh_rep_fyl):
        gencpy(sh_dir,workdir1,sh_rep_fyl)

    # Change to working directory
    os.chdir(workdir1)

    # Edit shell script files if needed (edit_sh_fyle = 1)
    if edit_sh_fyle:
        if inp_type == 'melts':
            edit_sh_files(workdir1,cont_run,biomass,inp_type,poly_conffile\
                          ,0,0,poly_topedit,'None','None',sh_pp_fyle\
                          ,sh_md_fyle,ff_dir,'None',0)
        elif inp_type == 'solvents':
            edit_sh_files(workdir1,cont_run,biomass,inp_type,poly_conffile,\
                          n_orgsolv,0,poly_topedit,o_sol_typ,'None',\
                          sh_pp_fyle,sh_md_fyle,ff_dir,sol_cfg,box_dim)
        elif inp_type == 'cosolvents':
            edit_sh_files(workdir1,cont_run,biomass,inp_type,poly_conffile,\
                          n_orgsolv,nwater,poly_topedit,o_sol_typ,wat_type,\
                          sh_pp_fyle,sh_md_fyle,ff_dir,sol_cfg,box_dim)

    # Write end of loop and return directory handle to main directory
    print( "End of run number: ", run_arr[casenum])
    os.chdir(main_dir)# current dir
