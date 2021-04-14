# Auxiliary data for generic plot
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import fileinput
import subprocess
import numpy as np
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

# Make labels and graphics


# Supporting files for input_gen.py
# Version: Mar-15-2021
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import fileinput
import subprocess
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

def ret_ana_dir(scr_dir,inp_type,biomass,disperse,casenum,subdir,\
                 solv_type='None'):
    # Check all directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        raise RuntimeError(head_dir + " does not exist!")

    fig_dir  = head_dir + '/results_figs'
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
        
    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        raise RuntimeError(poly_dir + " does not exist!")
        
    if inp_type == 'melts':
        poly_dir = poly_dir + '/' + disperse
        if not os.path.isdir(poly_dir):
            raise RuntimeError(poly_dir + " does not exist!")

    rundir = poly_dir + '/run_' + str(casenum)
    if not os.path.isdir(rundir):
        raise RuntimeError(rundir + " does not exist!")

    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, " does not exist!")

    anadir = workdir1 + '/' + subdir
    if not os.path.isdir(workdir1):
        raise RuntimeError(anadir, " does not exist!")

    return workdir1, anadir, fig_dir
#------------------------------------------------------------------        

# Change width of seaborn barplot
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # Change the bar width
        patch.set_width(new_value)

        # Recenter the bar
        patch.set_x(patch.get_x() + diff * .5)
#------------------------------------------------------------------
