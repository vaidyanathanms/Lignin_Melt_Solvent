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
o_sol_typ = 'EOH'  # relevant only for solvents/cosolvents
num_cases = 3      # for solvents, this will be == len(norg_solv)
nchains   = 1     # number of polymer chains
norg_solv = [1000] # number of organic solvents
nwater    = [4000] # number of water molecules (for cosolvents)
#------------------------------------------------------------------

# Directory Paths
main_dir  = os.getcwd()
gmx_dir   = '../src_gmx'
top_dir   = gmx_dir + '/solv_files/topol'
conf_dir  = gmx_dir + '/solv_files/initguess'
prm_dir   = gmx_dir + '/solv_files/prm_itp'
mdp_dir   = gmx_dir + '/' + 'mdp_files'
sh_dir    = gmx_dir + '/' + 'sh_files'
scr_dir   = '/gpfs/alpine/bip189/scratch/vaidyams/'
if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
scr_dir  = scr_dir + '/lignin'
if not os.path.isdir(scr_dir):
    os.mkdir(scr_dir)
#------------------------------------------------------------------

# Required GMX/sh Files
mdp_fyles  = ['minim_pyinp.mdp','nvt_pyino.mdp',\
              'npt_berendsen_pyinp.mdp','npt_main_pyinp.mdp']
sh_md_fyle = 'run_md_pyinp.inp'
sh_pp_fyle = 'run_preprocess_pyinp.sh'
#------------------------------------------------------------------

#Main Code
for casenum in range(num_cases):
    
    print( "Run number: ", casenum + 1)
      
    # Make directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        os.mkdir(head_dir)

    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        os.mkdir(poly_dir)

    if inp_type == 'melts':
        poly_dir = poly_dir + '/' + disperse
        if not os.path.isdir(poly_dir):
            os.mkdir(poly_dir)
    
    workdir1 = poly_dir + '/run_' + str(casenum+1)
    if not os.path.isdir(workdir1):
        print('ERR: Create path and input top/pdb files first')
        continue

    # Set thermostat variables (change if there is a temp sweep)
    Tetau_nvt,Tetau_berend,Tetau_parrah,Prtau_berend,Prtau_parrah,\
        ref_temp,ref_pres = couple_coeff(inp_type,coeff_fyle)

    # Copy and edit mdp files
    check_cpy_mdp_files(mdp_dir,workdir1,mdp_fyles,Tetau_nvt,\
                        Tetau_berend,Tetau_parrah,Prtau_berend,\
                        Prtau_parrah,ref_temp,ref_pres,tc_grps,tc_type)
    cont_run,edit_sh_fyle = cpy_edit_sh_files(sh_dir,workdir1,\
                                              sh_pp_fyle,sh_md_fyle)

    # Change to working directory
    os.chdir(workdir1)

    # Check for pdb/psf/top files
    check_inp_files(workdir1)

    # Copy/Edit top/conf/prm files for solvents/cosolvents
    ff_dir = create_ff_dir(workdir1)
    if inp_type == 'solvents' or inp_type == 'cosolvents':
        if cont_run = 0:
            top_fyl,prm_fyl=cpy_solv_files(top_dir,gro_dir,prm_dir,\
                                           ff_dir,inp_type,o_sol_typ,\
                                           wat_typ)

            edit_main_top_file(ff_dir,top_fyl,prm_fyl)

    # Edit shell script files

    
        
