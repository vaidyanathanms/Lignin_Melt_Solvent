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
def check_arr_dim(ot_len,on_len,wt_len,wn_len,no_len,nw_len,tlen\
                  ,hi_tlen,blen):
    if ot_len != on_len or ot_len != wt_len or ot_len != wn_len or \
       ot_len != no_len or ot_len != nw_len or ot_len != tlen \
                 or ot_len != hi_tlen or ot_len != blen:
        raise RuntimeError('Unequal input lengths')
#------------------------------------------------------------------

# Assign input values
def assign_vals(otyp_arr,oname_arr,wtyp_arr,wname_arr,norg_arr,\
                nwat_arr,box_arr,temp_arr,hi_t_arr,inp_type,ival):

    if inp_type == 'melts':
        o_sol_typ = 'None' 
        solv_name = 'None'
        wat_type  = 'None'
        wat_name  = 'None'
        n_orgsolv = 0
        nwater = 0
        ref_temp = temp_arr[ival]
        hi_ref_t = hi_t_arr[ival]
        box_dim  = box_arr[ival]
    elif inp_type == 'cosolvents':
        o_sol_typ = otyp_arr[ival]
        solv_name = oname_arr[ival]
        wat_type  = wtyp_arr[ival]
        wat_name  = wname_arr[ival]
        n_orgsolv = norg_arr[ival]
        nwater    = nwat_arr[ival]
        ref_temp  = temp_arr[ival]
        hi_ref_t = hi_t_arr[ival]
        box_dim   = box_arr[ival]
    elif inp_type == 'solvents':
        o_sol_typ = otyp_arr[ival]
        solv_name = oname_arr[ival]
        wat_type  = 'None'
        wat_name  = 'None'
        n_orgsolv = norg_arr[ival]
        nwater    = 0
        ref_temp  = temp_arr[ival]
        hi_ref_t  = hi_t_arr[ival]
        box_dim   = box_arr[ival]

    return o_sol_typ,solv_name,wat_type,wat_name,n_orgsolv,\
        nwater,ref_temp,hi_ref_t,box_dim
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
        
# Set default thermostat coefficients
def couple_coeff(inp_type,coeff_fyle = 'None'):
    # default. change if needed

    ref_pres = 1
    tau_temp_nvt     = 0.1
    tau_temp_nvthigh = 0.2
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
                if words[0] == 'TauTemp_NVT':
                    tau_temp_nvt  = float(words[1])
                elif words[0] == 'TauTempHigh_NVT':
                    tau_temp_nvthigh  = float(words[1])
                elif words[0] == 'TauTemp_Berendsen':
                    tau_temp_berend = float(words[1])
                elif words[0] == 'TauPres_Berendsen':
                    tau_pres_berend = float(words[1])
                elif words[0] == 'TauTemp_Parrinello':
                    tau_temp_parrah = float(words[1])
                elif words[0] == 'TauPres_Parrinello':
                    tau_pres_parrah = float(words[1])
                elif words[0] == 'Ref_Pres':
                    ref_pres  = float(words[1])
                elif words[0] == 'Melt_Topfile':
                    melt_topname = words[1]
                else:
                    raise RuntimeError("Unknown keyword: "+ words[0] \
                                       + " in " + str(coeff_fyle))
    return tau_temp_nvt,tau_temp_nvthigh,tau_temp_berend,\
        tau_temp_parrah,tau_pres_berend,tau_pres_parrah,ref_pres,\
        melt_topname
#------------------------------------------------------------------

#Check pdb/psf/top files for the melt (or polymer)
def check_inp_files(dum_inpdir,top_name):
    # check structure files (*.pdb/.gro)
    if glob.glob(dum_inpdir+'/*.pdb') == [] and \
       glob.glob(dum_inpdir+'/*.gro') == []:
        raise RuntimeError("No polymer pdb/gro files found")
    elif len(glob.glob(dum_inpdir+'/*.gro')) == 1:
        conf_fname = glob.glob(dum_inpdir+'/*.gro')[0]
    elif len(glob.glob(dum_inpdir+'/*.pdb')) == 1:
        conf_fname = glob.glob(dum_inpdir+'/*.pdb')[0]
    elif len(glob.glob(dum_inpdir+'/*.gro')) > 1:
        print('More than one config file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.gro')
        conf_fname = max(fnames, key=os.path.getctime)

    # check topology files
    if glob.glob(dum_inpdir+'/*.top') == []:
        raise RuntimeError("No polymer topology files found")
    if top_name != 'None':
        if not os.path.exists(dum_inpdir + '/' + top_name):
            raise RuntimeError("Specified poly top file not found")
        else:
            topol_fname = top_name
    elif len(glob.glob(dum_inpdir+'/*.top')) == 1 and \
         top_name == 'None':
        topol_fname = glob.glob(dum_inpdir+'/*.top')[0]
    else:
        print('More than one topology file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.top')
        topol_fname = max(fnames, key=os.path.getctime)
    
    return conf_fname, topol_fname

#------------------------------------------------------------------

# Create index groups
def create_indxgrps(destdir,inp_type,npoly_res,solv_name,wat_name):
    tcgrp_fname = 'tcgrp_inp.txt'
    if inp_type == 'melts':
        tc_str = 'System'
    elif inp_type == 'solvents':
        tc_str = 'resnr 1 to '+ str(npoly_res)
        tc_str = tc_str + '; ' + 'resname ' + solv_name + '\n'
    elif inp_type == 'cosolvents':
        tc_str = 'resnr 1 to '+ str(npoly_res)
        tc_str = tc_str + '; ' + 'resname ' + solv_name
        tc_str = tc_str + '; ' + 'resname ' + \
                 wat_name.split('_')[0] + '\n'
    
    with open(destdir+'/'+tcgrp_fname,'w') as fw:
        fw.write(tc_str)
    return tcgrp_fname
#------------------------------------------------------------------

# Create temperature coupling groups
def create_tcgrps(inp_type,npoly_res,solv_name,wat_name):
    if inp_type == 'melts':
        tc_grp = 'System'
        tc_typ = 'Single'
    elif inp_type == 'solvents':
        tc_grp = 'resnr_1_to_'+ str(npoly_res)
        tc_grp = tc_grp + '  ' + 'resname_' + solv_name
        tc_typ = 'Multi '
    elif inp_type == 'cosolvents':
        tc_grp = 'resnr_1_to_'+ str(npoly_res)
        tc_grp = tc_grp + '  ' + 'resname_' + solv_name
        tc_grp = tc_grp + '  ' + 'resname_' + wat_name
        tc_typ = 'Multi '
    return tc_grp, tc_typ
#------------------------------------------------------------------

# Check for mdp files and copy/edit if not present
def check_cpy_mdp_files(srcdir,destdir,mdp_fyles,inp_type,Tetau_nvt\
                        ,Tetau_highnvt,Tetau_berend,Tetau_parrah,\
                        Prtau_berend,Prtau_parrah,ref_temp,hi_ref_t,\
                        ref_pres,tc_grpdata,tc_grptype,headdir,coeff_fyle):

    if coeff_fyle != 'None': #gmx inp file from sys.argv
        print("copying", coeff_fyle)
        gencpy(headdir,destdir,coeff_fyle)

    # Only temperatures need to have multiple coupling groups
    str_Tetau_nvt   = genstr(inp_type,Tetau_nvt)
    str_Tetau_hinvt = genstr(inp_type,Tetau_highnvt)
    str_Tetau_ber   = genstr(inp_type,Tetau_berend)
    str_Tetau_par   = genstr(inp_type,Tetau_parrah)
    str_temp        = genstr(inp_type,ref_temp)
    str_hi_temp     = genstr(inp_type,hi_ref_t)

    for mdp_fname in mdp_fyles:
        if not os.path.exists(srcdir + '/' + mdp_fname):
            raise RuntimeError(mdp_fname," not found in ",srcdir)
        gencpy(srcdir,destdir,mdp_fname)
        py_fname = mdp_fname
        rev_fname = mdp_fname.replace('_pyinp','')
        fr  = open(destdir + '/' + py_fname,'r')
        fw  = open(destdir + '/' + rev_fname,'w')
        fid = fr.read().replace("py_tcgrps",tc_grpdata).\
              replace("py_grptype",tc_grptype).\
              replace("py_Temptau_vr",str_Tetau_nvt).\
              replace("py_HighTemptau_vr",str_Tetau_hinvt).\
              replace("py_Temptau_Berend", str_Tetau_ber).\
              replace("py_Temptau_ParRah",str_Tetau_par).\
              replace("py_Prestau_Berend",str(Prtau_berend)).\
              replace("py_Prestau_ParRah",str(Prtau_parrah)).\
              replace("py_ref_t",str_temp).\
              replace("py_Highref_t",str_hi_temp).\
              replace("py_ref_p",str(ref_pres))
        fw.write(fid)
        fw.close()
        fr.close()
#------------------------------------------------------------------

# Generate strings
def genstr(inp_type,inp_vals):
    if inp_type == 'melts':
        out_str = str(inp_vals)
    elif inp_type == 'solvents':
        out_str = str(inp_vals) + '  ' + str(inp_vals)
    elif inp_type == 'cosolvents':
        out_str = str(inp_vals) + '  ' + str(inp_vals) + '  ' + \
                  str(inp_vals)
    return out_str
#------------------------------------------------------------------

# Copy shell script files
def cpy_sh_files(srcdir,destdir,sh_pp_fyle,sh_md_fyle):

    # Check for tpr files to provide run conditions
    if glob.glob(destdir+'/*.tpr') == []:
        print('No tpr files found. Beginning new runs..')
        continue_run = 0
    else:
        print('tpr files found. Continuing runs..')
        continue_run = 1

    if continue_run == 0: 
        if not os.path.exists(srcdir+'/' + sh_pp_fyle):
            raise RuntimeError(sh_pp_fyle+" not found in "+\
                               srcdir)
        gencpy(srcdir,destdir,sh_pp_fyle)

        if not os.path.exists(srcdir+'/' + sh_md_fyle):
            raise RuntimeError(sh_md_fyle+" not found in "+\
                               srcdir)
        gencpy(srcdir,destdir,sh_md_fyle)
        

    if continue_run == 1: #at least one tpr file present
        sh_mdrun = sh_md_fyle.replace('_pyinp','')
        print("Recopying ", sh_md_fyle)
        gencpy(srcdir,destdir,sh_md_fyle)

    return continue_run
#------------------------------------------------------------------

# Create forcefield directory
def create_ff_dir(destdir):
    ff_dir = destdir + '/forcefields.ff'
    if not os.path.isdir(ff_dir):
        os.mkdir(ff_dir)
    return ff_dir
#------------------------------------------------------------------

# Check and copy solvent topology/initguess/itp files
def cpy_solv_files(top_dir,conf_dir,itp_dir,ff_dir,inp_type,\
                   osoltype,wat_typ='tip3p'):
    
    top_fyl = []; itp_fyl = []; conf_fyl = [];
    if inp_type == 'solvents' or inp_type == 'cosolvents':

        # Check topology
        src_dir = top_dir; o_file = osoltype + '.top'
        if not os.path.exists(src_dir + '/' + o_file):
            raise RuntimeError(o_file + " not found in "+src_dir)
        gencpy(src_dir,ff_dir,o_file)
        top_fyl.append(o_file)

        # Check itp
        src_dir = itp_dir; o_file = osoltype + '.itp'
        if not os.path.exists(src_dir + '/' + o_file):
            raise RuntimeError(o_file+" not found in " + src_dir)
        gencpy(src_dir,ff_dir,o_file)
        itp_fyl.append(o_file)

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
            raise RuntimeError(osoltype+" conf not found in "+src_dir)

    if inp_type == 'cosolvents': #with second solvent; def: water
        
        # check top/itp/conf for cosolvent
        src_dir = top_dir; w_file = wat_typ + '.top'
        if not os.path.exists(src_dir + '/' + w_file):
            raise RuntimeError(w_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,w_file)
        top_fyl.append(w_file)

        src_dir = itp_dir; w_file = wat_typ + '.itp'
        if not os.path.exists(src_dir + '/' + w_file):
            raise RuntimeError(w_file + "not found in " + src_dir)
        gencpy(src_dir,ff_dir,w_file)
        itp_fyl.append(w_file)

        src_dir = conf_dir
        w_file1 = wat_typ + '.gro'; w_file2 = wat_typ + '.pdb'; 
        if os.path.exists(src_dir + '/' + w_file1):
            gencpy(src_dir,ff_dir,w_file1)
            conf_fyl.append(w_file1)
        elif os.path.exists(src_dir + '/' + w_file2):
            gencpy(src_dir,ff_dir,w_file2)
            conf_fyl.append(w_file2)
        else:
            raise RuntimeError(wat_typ+" conf not found in "+src_dir)

    return top_fyl,itp_fyl,conf_fyl
#------------------------------------------------------------------

# Edit melt/polymer topology to add solvent/cosolvent details
def edit_main_top_file(main_topfyle,ff,top_arr,itp_arr,workdir):

    # check whether main_topfyle has the required data. If so no need
    # to edit and make the new ones. Slightly inefficient - but OK!
    req_flag = 0 # default not present
    for i in range(len(top_arr)):
        with open(main_topfyle) as fchk:
            for line in fchk:
                if top_arr[i] in line:
                    req_flag += 1 # found string

    for i in range(len(itp_arr)):
        with open(main_topfyle) as fchk:
            for line in fchk:
                if itp_arr[i] in line:
                    req_flag += 1 # found string

    sum_req = len(itp_arr) + len(top_arr)
    if req_flag == sum_req:
        print('Present topol file in the directory has all reqd data')
        return main_topfyle
    else:
        print('Adding reqd data to new topology file')

    # split strings for ease of reading in final file o/p
    ff_dir_spl = ff.split("/")
    ff_dir = "./" + ff_dir_spl[len(ff_dir_spl)-1]

    spl_str = main_topfyle.split("/")
    edited_file = workdir + "/solvated_" + spl_str[len(spl_str)-1]

    # create new topology file with required details
    # add topology before [ moleculetype ]
    inc_top = ''
    inc_pre = '#include '
    for i in range(len(top_arr)):
        inc_top = inc_top + inc_pre + "\"" + ff_dir + '/' + \
                  top_arr[i] + "\"" + '\n'
    inc_top = inc_top + '\n'

    # add parameters before [ system ]
    inc_itp = ''
    inc_pre = '#include '
    for i in range(len(itp_arr)):
        inc_itp = inc_itp + inc_pre + "\"" + ff_dir + '/' + \
                  itp_arr[i] + "\"" + '\n'
    inc_itp = inc_itp + '\n'

    with open(main_topfyle,"r") as fin:
        with open(edited_file,"w") as fout:
            for line in fin:
                if "[ moleculetype ]" in line:
                    line=line.replace(line,inc_itp+line)
                elif "[ system ]" in line:
                    line=line.replace(line,inc_top+line)
                fout.write(line)

    return edited_file
#------------------------------------------------------------------

# Edit pre-processing shell script files
def edit_pp_files(biomass,inp_type,polycfg,nsolv,nwater,\
                  topfyle,o_sol_type,wat_typ,pp_fyle,ffdir,\
                  solcfg,dim,indx_fyle):

    # use the last element (others correspond to full path)
    # underscored variable corresponds to split file name
    top_fyle = ret_file_str(topfyle)
    poly_cfg = ret_file_str(polycfg)
    ff_dir   = ret_file_str(ffdir)

    # job/box name
    jname = 'pp_'+ biomass
    if inp_type == 'solvents' or inp_type == 'cosolvents':
        jname = jname + '_' + o_sol_type #py_jobname
    box_conffyle = "boxedit_" + poly_cfg
    fin_conf = 'initconf.gro' # in gro format

    # solvate commands
    solv_js = 'jsrun -X 1 -n 1 -c 7 -a 1 -g 1 ' + \
              '--launch_distribution plane:1 ' + \
              '-b packed:7 gmx_mpi solvate'
    if inp_type == 'melts':
        solv_str1 = '# no solvation'
        solv_str2 = '# no cosolvation'
    elif inp_type == 'solvents':
        sol_cfg  = ret_file_str(solcfg[0])
        sol_cfg1 = ff_dir + '/' + sol_cfg
        solv_str1 = solv_js + ' -cp ' + poly_cfg + ' -cs ' + \
                    sol_cfg1 + ' -p ' + top_fyle + ' -o ' + \
                    fin_conf + ' -maxsol ' + str(nsolv)\
                    + ' -box ' + str(dim) + ' ' + str(dim) + ' '\
                    + str(dim) + '\n'
        solv_str2 = '# no cosolvation'
    elif inp_type == 'cosolvents':
        sol_cfg  = ret_file_str(solcfg[0])
        sol_cfg1 = ff_dir + '/' + sol_cfg
        sol_cfg  = ret_file_str(solcfg[1])
        sol_cfg2 = ff_dir + '/' + sol_cfg
        solv_str1 = solv_js + ' -cp ' + poly_cfg + ' -cs ' + \
                    sol_cfg1 + ' -p ' + top_fyle + ' -o ' + \
                    'solv_'+ poly_cfg + ' -maxsol ' + str(nsolv)\
                    + ' -box ' + str(dim) + ' ' + str(dim) + ' '\
                    + str(dim) + '\n'

        dim += 1 # arbitrary 1 nm
        solv_str2 = solv_js + ' -cp ' +'solv_' +poly_cfg +' -cs '\
                    + sol_cfg2 + ' -p ' + top_fyle + ' -o ' + \
                    fin_conf +' -maxsol '+str(nwater)\
                    + ' -box ' + str(dim) + ' ' + str(dim) + ' '\
                    + str(dim) + '\n'

    # edit pp_fyle
    py_fname = pp_fyle
    rev_fname = py_fname.replace('_pyinp','')
    fr  = open(py_fname,'r')
    fw  = open(rev_fname,'w')
    fid = fr.read().replace("py_jobname",jname).\
          replace("py_meltconf",poly_cfg).\
          replace("py_boxmeltconf",box_conffyle).\
          replace("py_solvate_1",solv_str1).\
          replace("py_solvate_2", solv_str2).\
          replace("py_topol",top_fyle).\
          replace("py_indexfyle",indx_fyle).\
          replace("py_finconf",fin_conf)
    fw.write(fid)
    fr.close(); fw.close()

#------------------------------------------------------------------

# returns the last element file string after stripping "/"
def ret_file_str(inpstr):
    spl_str = inpstr.split("/")
    return spl_str[len(spl_str)-1]
#------------------------------------------------------------------

# Edit run_md.sh file all times
def edit_md_files(biomass,inp_type,polycfg,topfyle,o_sol_type\
                  ,md_fyle,indx_fyle):

    # use the last element (others correspond to full path)
    # underscored variable corresponds to split file name
    top_fyle = ret_file_str(topfyle)
    poly_cfg = ret_file_str(polycfg)

    # edit run_md_file ALWAYS
    # job/box name
    jname = 'md_'+ biomass
    if inp_type == 'solvents' or inp_type == 'cosolvents':
        jname = jname + '_' + o_sol_type #py_jobname

    fin_conf = 'initconf.gro' # in gro format

    # edit md_fyle
    py_fname = md_fyle
    rev_fname = py_fname.replace('_pyinp','')
    fr  = open(py_fname,'r')
    fw  = open(rev_fname,'w')
    fid = fr.read().replace("py_jobname",jname).\
          replace("py_indexfyle",indx_fyle).\
          replace("py_topol",top_fyle).\
          replace("py_finconf",fin_conf)
    fw.write(fid)
    fr.close(); fw.close()
#------------------------------------------------------------------
