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

# Check consistency of input dimensions
def check_arr_dim(ot_len,on_len,wt_len,wn_len):
    if ot_len != on_len or ot_len != wt_len or ot_len != wn_len:
        raise RuntimeError('Unequal input lengths')
#------------------------------------------------------------------

# Assign input values
def assign_vals(otyp_arr,oname_arr,wtyp_arr,wname_arr,inp_type,ival):

    if inp_type == 'melts':
        o_sol_typ = 'None' 
        solv_name = 'None'
        wat_type  = 'None'
        wat_name  = 'None'
    elif inp_type == 'cosolvents':
        o_sol_typ = otyp_arr[ival]
        solv_name = oname_arr[ival]
        wat_type  = wtyp_arr[ival]
        wat_name  = wname_arr[ival]
    elif inp_type == 'solvents':
        o_sol_typ = otyp_arr[ival]
        solv_name = oname_arr[ival]
        wat_type  = 'None'
        wat_name  = 'None'

    return o_sol_typ,solv_name,wat_type,wat_name
#------------------------------------------------------------------

# Set working directory
def set_working_dir(rundir,inp_type,solv_type = 'None'):
    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, "does not exist")
    
    return workdir1
#------------------------------------------------------------------
# Check and copy all required files: NEED TO WORK FOR MELTS
def check_cpy_ana_files(req_files,sh_fyle,inp_type,o_solv,\
                        ana_dir,workdir1):

    # Check all req_files and copy
    if inp_type == 'solvents' or inp_type == 'cosolvents':
        src_dir = ana_dir+'/'+o_solv
    elif inp_type == 'melts':
        src_dir = ana_dir+'/lig_melts'

    for fname in req_files:
        if not os.path.exists(src_dir+'/'+fname):
            raise RuntimeError(fname," not found in ",src_dir)
        else:
            gencpy(src_dir,workdir1,fname)

    # Check shell script file and copy
    if not os.path.exists(src_dir+'/'+sh_fyle):
        raise RuntimeError(sh_fyle, " not found in ", src_dir)
    else:
        gencpy(src_dir,workdir1,sh_fyle)
#------------------------------------------------------------------

# Retrieve trr/tpr files
def retrieve_traj_files(workdir):
    os.chdir(workdir)
    tpr_list = glob.glob('/*.tpr')
    if tpr_list == []:
        raise RuntimeError('No tpr files found in ', workdir)
    
    trr_list = glob.glob('/*.trr')
    if trr_list == []:
        raise RuntimeError('No trr files found in ', workdir)

    print('Retrieving latest tpr/trr files')
    tpr_fyle = max(tpr_list,key=os.path.getctime)
    trr_fyle = max(tpr_list,key=os.path.getctime)
    
    return tpr_fyle, trr_fyle
#------------------------------------------------------------------
# TO BE COMPLETED-to read from external file what all needs to be
# analyzed 
def parse_and_write(anainp,run_sh,workdir):

    os.chdir(workdir)
    if not os.path.exists(anainp):
        raise RuntimeError(anainp, " not found in ", workdir)
    
    with open(anainp) as fin:
        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('#'):
                continue
            words = line.split()
#------------------------------------------------------------------
# Edit runana file
def edit_anash_file(biomass,inp_type,o_sol_type,sh_anafyl):

    # edit run_ana file
    # job/box name
    jname = 'ana_'+ biomass
    if inp_type == 'solvents' or inp_type == 'cosolvents':
        jname = jname + '_' + o_sol_type #py_jobname

    fin_conf = 'initconf.gro' # in gro format

    # edit md_fyle
    py_fname = sh_anafyl
    rev_fname = py_fname.replace('_pyinp','')
    fr  = open(py_fname,'r')
    fw  = open(rev_fname,'w')
    fid = fr.read().replace("py_jobname",jname)
    fw.write(fid)
    fr.close(); fw.close()
#------------------------------------------------------------------
