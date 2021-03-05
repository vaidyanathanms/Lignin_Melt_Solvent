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
def couple_coeff(inp_type,coeff_fyle = 'None'):
    # default. change if needed
    tau_temp_nvt     = 0.1
    tau_temp_berend  = 0.1
    tau_temp_parrah  = 0.2
    tau_pres_berend  = 0.5
    ref_temp = 300
    ref_pres = 1
    if inp_type == 'melts':
        tau_pres_parrah  = 5.0
    else:
        tau_pres_parrah  = 8.0
    if coeff_fyle != 'None':
        with open(coeff_fyle) as farg:
            for line in farg:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    continue
                words = line.split()
                if words[0] = 'TauTemp_NVT':
                    tau_temp_nvt  = 0.1
                elif words[0] == 'TauTemp_Berendsen':
                    tau_temp_berend = float(words[1])
                elif words[0] == 'TauTemp_Parrinello':
                    tau_temp_parrah = float(words[1])
                elif words[0] == 'TauPres_Berendsen':
                    tau_pres_berend = float(words[1])
                elif words[0] == 'TauPres_Parrinello':
                    tau_pres_parrah = float(words[1])
                elif words[0] == 'Ref_Temp':
                    ref_temp  = float(words[1])
                elif words[0] == 'Ref_Pres':
                    ref_pres  = float(words[1])
                else:
                    raise RuntimeError("Unknown keyword: "+ words[0] \
                                       + "in" + str(coeff_fyle))
    return tau_temp_nvt, tau_temp_berend, tau_temp_parrah, \
        tau_pres_berend, tau_pres_parrah, ref_temp, ref_pres
#------------------------------------------------------------------

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
def check_cpy_mdp_files(srcdir,destdir,mdp_fyles,Tetau_nvt,\
                        Tetau_berend,Tetau_parrah,Prtau_berend,\
                        Prtau_parrah,ref_temp,ref_pres,tc_grpdata,\
                        tc_grptype):
     if not os.path.exists(mdp_fyles[fylcnt]):
         gencpy(srcdir,desdir,mdp_fyles[fylcnt])
         py_fname = mdp_fyles[fylcnt]
         rev_fname = py_fname.replace('_pyinp','')
         print(rev_fname)
         fr  = open(py_fname,'r')
         fw  = open(rev_fname,'w')

         fid = fr.read().replace("py_tcgrps",tc_grpdata).\
               replace("py_grptype",tc_grptype).\
               replace("py_Temptau_vr",str(Tetau_nvt)).\
               replace("py_Temptau_Berend", str(Tetau_berend)).\
               replace("py_Temptau_ParRah",str(Tetau_parrah)).\
               replace("py_Prestau_Berend",str(Prtau_berend)).\
               replace("py_Prestau_ParRah", str(Prtau_parrah)).\
               replace("py_ref_t",str(ref_temp)).\
               replace("py_ref_p",str(ref_pres))
         fw.write(fid)
         fw.close()
         fr.close()

#------------------------------------------------------------------
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
            if not os.path.exists(sh_fyles[fylcnt]):
                raise RuntimeError(sh_fyles[fylcnt] + "not found in" \
                                   + srcdir)
            gencpy(srcdir,desdir,sh_fyles[fylcnt])
            





#Main Code
