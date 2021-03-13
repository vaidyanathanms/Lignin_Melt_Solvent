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
    ref_temp = 300
    ref_pres = 1
    tau_temp_nvt     = 0.1
    tau_temp_berend  = 0.1
    tau_temp_parrah  = 0.2
    tau_pres_berend  = 0.5
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
                    tau_temp_nvt  = float(words[1])
                elif words[0] == 'TauTemp_Berendsen':
                    tau_temp_berend = float(words[1])
                elif words[0] == 'TauPres_Berendsen':
                    tau_pres_berend = float(words[1])
                elif words[0] == 'TauTemp_Parrinello':
                    tau_temp_parrah = float(words[1])
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
def check_inp_files(dum_inpdir,top_name = 'None'):
    # check structure files (.pdb/.gro)
    if glob.glob(dum_inpdir+'/*.pdb') == [] and \
       glob.glob(dum_inpdir+'/*.gro') == []:
        raise RuntimeError("No polymer pdb/gro files found")
    if glob.glob(dum_inpdir+'/*.top') == []:
        raise RuntimeError("No polymer topology files found")
    if top_name != 'None':
        if not os.path.exists(dum_inpdir + '/' + top_name):
            raise RuntimeError("Specified poly top file not found")
    else len(glob.glob(dum_inpdir+'/*.top')) > 1 and \
       top_name == 'None':
        print('More than one topology file found. Using latest')
           
#------------------------------------------------------------------

# Check for mdp files and copy/edit if not present
def check_cpy_mdp_files(srcdir,destdir,mdp_fyles,Tetau_nvt,\
                        Tetau_berend,Tetau_parrah,Prtau_berend,\
                        Prtau_parrah,ref_temp,ref_pres,tc_grpdata,\
                        tc_grptype,coeff_fyle='None'):

    if coeff_fyle != 'None':
        gencpy(srcdir,destdir,coeff_fyle)

    if not os.path.exists(mdp_fyles[fylcnt]):
        gencpy(srcdir,destdir,mdp_fyles[fylcnt])
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

# Copy and edit shell script files
def cpy_edit_sh_files(srcdir,destdir,sh_pp_fyle,sh_md_fyle):
    edit_sh_fyle = 0
    # Check for tpr files to provide run conditions
    if glob.glob(destdir+'/*.tpr') == []:
        print('No tpr files found. Beginning new runs..')
        continue_run = 0
        edit_sh_fyle = 1
    else:
        print('tpr files found. Continuing runs..')
        continue_run = 1

    if rerun == 0: 
        if not os.path.exists(srcdir+'/' + sh_pp_fyle):
            raise RuntimeError(sh_pp_fyle + "not found in" +\
                               srcdir)
        gencpy(srcdir,desdir,sh_pp_fyle)

        if not os.path.exists(srcdir+'/' + sh_md_fyle):
            raise RuntimeError(sh_md_fyle + "not found in" +\
                               srcdir)
        gencpy(srcdir,desdir,sh_md_fyle)
        

    if rerun == 1: #at least one tpr file present
        sh_mdrun = sh_md_fyle.replace('_pyinp','')
        if not os.path.exists(destdir+'/' + sh_mdrun):
            print("WARNING: ", sh_mdrun, "not found in", destdir)
            print("Recopying ", sh_md_fyle)
            edit_sh_fyle = 1
            gencpy(srcdir,desdir,sh_md_fyle)
            
    return continue_run, edit_sh_fyle
#------------------------------------------------------------------

# Create forcefield directory
def create_ff_dir(destdir):
    ff_dir = destdir + '/forcefields.ff'
    if not os.path.isdir(ff_dir):
        os.mkdir(ff_dir)
    return ff_dir
#------------------------------------------------------------------

# Check and copy solvent topology/initguess/prm files
def cpy_solv_files(top_dir,conf_dir,prm_dir,destdir,inp_type,\
                   osoltype,wat_typ='tip3p'):
    
    top_fyl = []; prm_fyl = []
    if inp_type == 'cosolvents':

        # Check topology
        src_dir = top_dir
        ext = '.top'
        o_file = osoltype + ext
        if not os.path.exists(src_dir + '/' + o_file):
            raise RuntimeError(o_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,o_file)
        top_fyl.append(o_file)

        w_file = wat_typ + ext
        if not os.path.exists(src_dir + '/' + w_file):
            raise RuntimeError(w_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,w_file)
        top_fyl.append(w_file)

        # Check prm
        src_dir = prm_dir
        ext = '.prm'
        o_file = osoltype + ext
        if not os.path.exists(src_dir + '/' + o_file):
            raise RuntimeError(o_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,o_file)
        prm_fyl.append(o_file)

        w_file = wat_typ + ext
        if not os.path.exists(src_dir + '/' + w_file):
            raise RuntimeError(w_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,w_file)
        prm_fyl.append(w_file)

        # Check gro/pdb conf file
        src_dir = conf_dir
        ext1 = '.gro'; ext2 = '.pdb'
        o_file1 = osoltype + ext1; o_file2 = osoltype + ext2; 
        if os.path.exists(src_dir + '/' + o_file1):
            gencpy(src_dir,ff_dir,o_file1)
        elif os.path.exists(src_dir + '/' + o_file2):
            gencpy(src_dir,ff_dir,o_file2)
        else:
            raise RuntimeError(osoltype+"conf not found in "+src_dir)

        src_dir = conf_dir
        ext1 = '.gro'; ext2 = '.pdb'
        w_file1 = wat_typ + ext1; w_file2 = wat_typ + ext2; 
        if os.path.exists(src_dir + '/' + w_file1):
            gencpy(src_dir,ff_dir,w_file1)
        elif os.path.exists(src_dir + '/' + w_file2):
            gencpy(src_dir,ff_dir,w_file2)
            w_conf_fyl = '.pdb'
        else:
            raise RuntimeError(wat_typ+"conf not found in "+src_dir)

    return top_fyl,prm_fyl

def edit_main_top_file(main_topfyle,ff_dir,top_fyl,prm_fyl):
    
