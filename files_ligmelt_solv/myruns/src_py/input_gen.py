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

    # Set thermostat variables (change if there is a temp sweep)
    Tetau_nvt,Tetau_berend,Tetau_parrah,Prtau_berend,Prtau_parrah,\
        ref_temp, ref_pres = couple_coeff(inp_type,coeff_fyle)

    # Change to working directory
    os.chdir(workdir1)

    # Check for pdb/psf/top files
    check_inp_files(workdir1)

    # Copy and edit mdp files
    check_cpy_mdp_files(mdp_dir,workdir1,mdp_fyles,Tetau_nvt,\
                        Tetau_berend,Tetau_parrah,Prtau_berend,\
                        Prtau_parrah,ref_temp,ref_pres,tc_grps,tc_type)
    rerun = cpy_edit_sh_files(sh_dir,workdir1,sh_fyles)

    
        
