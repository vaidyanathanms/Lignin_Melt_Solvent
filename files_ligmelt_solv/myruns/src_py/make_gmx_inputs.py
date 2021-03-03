# Supporting files for input_gen.py
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
#------------------------------------------------------------------

# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
#------------------------------------------------------------------

# Set default thermostat coefficients
def couple_coeff(coeff_fyle = 'None'):
    temp_velresc = 0.1
    temp_parrah  = 0.3
    pres_parrah  = 8.0
    if coeff_fyle = 'None':
        

#Check pdb/psf/top files
def check_inp_files(dum_inpdir):
    # check structure files (.pdb/.gro)
    if glob.glob(dum_inpdir+'/*.pdb') == [] and \
       glob.glob(dum_inpdir+'/*.gro') == []:
        raise RuntimeError("No structure (pdb/gro) files found")
    if glob.glob(dum_inpdir+'/*.top') == []:
        raise RuntimeError("No topology (top) files found")
#------------------------------------------------------------------

# Check for mdp files and copy/edit if not present
def check_cpy_mdp_files(srcdir,destdir,mdp_fyles):


# Copy and edit mdp files
def cpy_edit_sh_files(srcdir,destdir,sh_fyles):
    # Check for tpr files to provide run conditions
    if glob.glob(destdir+'/*.tpr') == []:
        print('No tpr files found. Beginning new runs..')
        rerun = 0
    else:
        print('tpr files found. Continuing runs..')
        rerun = 1

    if rerun == 1: #tpr files present
        for fylcnt in range(len(sh_fyles)):
            if not os.path.exists(mdp_fyles[fylcnt]):
                raise RuntimeError(mdp_fyles[fylcnt] + "not found in" \
                                   + srcdir)
            gencpy(srcdir,desdir,mdp_fyles[fylcnt])
            


#Directory Paths


#Required Files


#Main Code
