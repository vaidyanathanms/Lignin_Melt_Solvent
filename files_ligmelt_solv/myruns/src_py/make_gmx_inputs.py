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
import fileinput
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

# Set working directory
def set_working_dir(rundir,inp_type,solv_type = 'None'):
    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        workdir1 = h_workdir + '/' + solv_type + '_water'
    elif:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, "does not exist")
    
    return workdir1
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
    melt_topname     = 'None'
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
                elif words[0] == 'Melt_Topfile':
                    melt_topname = words[1]
                else:
                    raise RuntimeError("Unknown keyword: "+ words[0] \
                                       + "in" + str(coeff_fyle))
    return tau_temp_nvt,tau_temp_berend, tau_temp_parrah, \
        tau_pres_berend,tau_pres_parrah,ref_temp,ref_pres, \
        melt_topname
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

# Copy shell script files
def cpy_sh_files(srcdir,destdir,sh_pp_fyle,sh_md_fyle):
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

#Check pdb/psf/top files for the melt (or polymer)
def check_inp_files(dum_inpdir,top_name):
    # check structure files (.pdb/.gro)
    if glob.glob(dum_inpdir+'/*.pdb') == [] and \
       glob.glob(dum_inpdir+'/*.gro') == []:
        raise RuntimeError("No polymer pdb/gro files found")
    elif len(glob.glob(dum_inpdir+'/*.gro')) == 1:
        conf_name = glob.glob(dum_inpdir+'/*.gro')
    elif len(glob.glob(dum_inpdir+'/*.pdb')) == 1:
        conf_name = glob.glob(dum_inpdir+'/*.pdb')
    elif len(glob.glob(dum_inpdir+'/*.gro')) > 1:
        print('More than one config file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.gro')
        conf_fname = max(fnames, key=os.path.getctime)
    elif len(glob.glob(dum_inpdir+'/*.gro')) > 1:
        print('More than one config file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.pdb')
        conf_fname = max(fnames, key=os.path.getctime)
        

    # check topology files
    if glob.glob(dum_inpdir+'/*.top') == []:
        raise RuntimeError("No polymer topology files found")
    if top_name != 'None':
        if not os.path.exists(dum_inpdir + '/' + top_name):
            raise RuntimeError("Specified poly top file not found")
        else:
            topol_name = top_name
    elif len(glob.glob(dum_inpdir+'/*.top')) == 1 and \
         top_name == 'None':
        topol_name = glob.glob(dum_inpdir+'/*.top')
    else:
        print('More than one topology file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.top')
        topol_fname = max(fnames, key=os.path.getctime)
    
    return conf_fname, topol_fname

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
    
    top_fyl = []; prm_fyl = []; conf_fyl = [];
    if inp_type == 'solvents' or inp_type == 'cosolvents':

        # Check topology
        src_dir = top_dir; o_file = osoltype + '.top'
        if not os.path.exists(src_dir + '/' + o_file):
            raise RuntimeError(o_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,o_file)
        top_fyl.append(o_file)

        # Check prm
        src_dir = prm_dir; o_file = osoltype + '.prm'
        if not os.path.exists(src_dir + '/' + o_file):
            raise RuntimeError(o_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,o_file)
        prm_fyl.append(o_file)

        # Check gro/pdb conf file
        src_dir = conf_dir
        o_file1 = osoltype + '.gro'; o_file2 = osoltype + '.pdb'
        if os.path.exists(src_dir + '/' + o_file1):
            gencpy(src_dir,ff_dir,o_file1)
            conf_fyl.append(o_file1)
        elif os.path.exists(src_dir + '/' + o_file2):
            gencpy(src_dir,ff_dir,o_file2)
            conf_fyl.append(o_file2)
        else:
            raise RuntimeError(osoltype+"conf not found in "+src_dir)

    if inp_type == 'cosolvents': #with second solvent; def: water
        
        # check top/prm/conf for cosolvent
        src_dir = top_dir; w_file = wat_typ + '.top'
        if not os.path.exists(src_dir + '/' + w_file):
            raise RuntimeError(w_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,w_file)
        top_fyl.append(w_file)

        src_dir = prm_dir; w_file = wat_typ + '.prm'
        if not os.path.exists(src_dir + '/' + w_file):
            raise RuntimeError(w_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,w_file)
        prm_fyl.append(w_file)

        src_dir = conf_dir
        w_file1 = wat_typ + '.gro'; w_file2 = wat_typ + '.pdb'; 
        if os.path.exists(src_dir + '/' + w_file1):
            gencpy(src_dir,ff_dir,w_file1)
            conf_fyl.append(w_file1)
        elif os.path.exists(src_dir + '/' + w_file2):
            gencpy(src_dir,ff_dir,w_file2)
            conf_fyl.append(w_file2)
        else:
            raise RuntimeError(wat_typ+"conf not found in "+src_dir)

    return top_fyl,prm_fyl,conf_fyl
#------------------------------------------------------------------

# Edit melt/polymer topology to add solvent/cosolvent details
def edit_main_top_file(workdir,main_topfyle,ff_dir,top_arr,prm_arr):
    # add topology before [ moleculetype ]
    inc_top = ''
    inc_pre = '#include '
    for i in len(range(top_arr)):
        inc_top = inc_top + inc_pre + "'" + ff_dir + '/' + \
                  top_fyl[i] + "' " + '/n'

    with open(workdir+'/'+main_topfyle) as fyl:
        for line in fyl.readlines():
            if "[ moleculetype ]" in line:
                line=line.replace(line,inc_top+line)

    # add parameters before [ system ]
    inc_prm = ''
    inc_pre = '#include '
    for i in len(range(prm_arr)):
        inc_prm = inc_top + inc_pre + "'" + ff_dir + '/' + \
                  top_fyl[i] + "' " + '/n'

    with open(workdir+'/'+main_topfyle) as fyl:
        for line in fyl.readlines():
            if "[ system ]" in line:
                line=line.replace(line,inc_top+line)
#------------------------------------------------------------------

# Edit shell script files
def edit_sh_files(workdir,cont_run,biomass,inp_type,poly_cfg,\
                  nsolv,nwater,topfyle,o_sol_type,wat_typ,pp_fyle,\
                  md_fyle,ff_dir,sol_cfg,dim):

    if cont_run == 0:
        # job/box name
        jname = 'pp_'+ biomass
        if inp_type == 'solvents' or inp_type == 'cosolvents':
            jname = jname + '_' + o_sol_type #py_jobname
        box_conffyle = "boxedit_" + conf_fyle
        
        # solvate commands
        solv_js = 'jsrun -X 1 -n 1 -c 7 -a 1 -g 1 ' + \
                  '--launch_distribution plane: 1 ' + \
                  '-b packed:7 gmx_mpi solvate'
        if inp_type == 'melts':
            solv_str1 = '# no solvation'
            solv_str2 = '# no cosolvation'
            fin_conf  = poly_cfg
        elif inp_type == 'solvents':
            sol_cfg1 = ff_dir + '/' + sol_cfg[0]
            solv_str1 = solv_js + ' -cp ' + poly_cfg + ' -cs ' + \
                        sol_cfg1 + ' -p ' + topfyle + ' -o ' + \
                        'solv_'+ poly_cfg + ' -maxsol ' + str(nsolv)\
                        + ' -box ' + str(dim) + ' ' + str(dim) + ' '\
                        + str(dim) + '\n'
            solv_str2 = '# no cosolvation'
            fin_conf  = 'solv_' + poly_cfg
        elif inp_type == 'cosolvents':
            sol_cfg1 = ff_dir + '/' + sol_cfg[0]
            sol_cfg2 = ff_dir + '/' + sol_cfg[1]
            solv_str1 = solv_js + ' -cp ' + poly_cfg + ' -cs ' + \
                        sol_cfg1 + ' -p ' + topfyle + ' -o ' + \
                        'solv_'+ poly_cfg + ' -maxsol ' + str(nsolv)\
                        + ' -box ' + str(dim) + ' ' + str(dim) + ' '\
                        + str(dim) + '\n'
            solv_str2 = solv_js + ' -cp ' +'solv_' +poly_cfg +' -cs '\
                        + sol_cfg2 + ' -p ' + topfyle + ' -o ' + \
                        'cosolv_'+ poly_cfg +' -maxsol '+str(nwater)\
                        + ' -box ' + str(dim) + ' ' + str(dim) + ' '\
                        + str(dim) + '\n'
            fin_conf  = 'cosolv_' + poly_cfg

        # edit pp_fyle
        py_fname = pp_fyle
        rev_fname = py_fname.replace('_pyinp','')
        fr  = open(py_fname,'r')
        fw  = open(rev_fname,'w')
        fid = fr.read().replace("py_meltconf",conf_fyle).\
              replace("py_boxmeltconf",box_conffyle).\
              replace("py_solvate_1",solv_str1).\
              replace("py_solvate_2", solv_str2).\
              replace("py_topol",topfyle).\
              replace("py_finconf",fin_conf)
        fr.close(); fw.close()

    # since edit_sh_fyle is true, the run_md_file is recopied
    # editing run_md_file
    # job/box name
    jname = 'md_'+ biomass
    if inp_type == 'solvents' or inp_type == 'cosolvents':
        jname = jname + '_' + o_sol_type #py_jobname
        box_conffyle = "boxedit_" + conf_fyle

    # edit pp_fyle
    py_fname = pp_fyle
    rev_fname = py_fname.replace('_pyinp','')
    fr  = open(py_fname,'r')
    fw  = open(rev_fname,'w')
    fid = fr.read().replace("py_topol",topfyle).\
          replace("py_finconf",fin_conf)
    fr.close(); fw.close()
#------------------------------------------------------------------
