# Generic plot
# Version: Mar-16-2021
#------------------------------------------------------------------

# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import sys
import re
import shutil
import glob
import math
import subprocess
from plt_aux import *
#------------------------------------------------------------------

# Color/line data; figure defaults
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 20}) # other fonts
#------------------------------------------------------------------

# Input data
inp_type  = 'cosolvents' # melts, solvents, cosolvents
disperse  = 'mono' # mono/poly; only for melts
biom_arr  = ['COMT'] # biomass type arr
otyp_arr  = ['EOH','GVL','THF']  # solvent arr for solvents/cosolvents
run_arr   = [4,5,6] # number of independent runs for a given biomass
#------------------------------------------------------------------

# Plot keys 0,1,2 (0-None,1-separate,2-together)
sasa    = 0 # plot SASA distribution and average
rdf_pol = 0 # plot polymer-solvent/water RDF 0 or 1
rdf_bO4 = 0 # plot bO4-solvent/water RDF
hbonds  = 0 # hbonds of solv/water, total hbonds, hbond dist
rad_gyr = 2 # avg rg and rg distribution
#------------------------------------------------------------------

# Directory paths
main_dir = os.getcwd() # current dir
scr_dir  = '/gpfs/alpine/bip189/scratch/vaidyams' # scratch dir
scr_dir  = scr_dir + '/lignin'
if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
#------------------------------------------------------------------

# Check basic directories and return head directory

# Plot avg SASA for all solvents
if sasa != 0:
    print("Analyzing SASA data")

    if sasa == 2: # plot all averages together
        consol_dir = scr_dir + '/' + inp_type + '/consolidated'
        if not os.path.isdir(consol_dir):
            os.mkdir(consol_dir)

        # Write data and then read to concatenate
        sasa_fyl = consol_dir + '/sasadata.dat'
        fc_ss = open(sasa_fyl,'w')
        fc_ss.write('%s\t%s\t%s\n' %('Biomass','Solvent','SASA'))

    fig1,ax1 = plt.subplots()
    ax1.set_xlabel(r'SASA (nm$^2$)')
    ax1.set_ylabel('Probability')
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        ysasa_avg = np.zeros(0)

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            yall = np.zeros(0)
            yavg = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  solv_type)
                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file
                if os.path.exists(anadir + '/sasa.xvg'):
                    fname  = anadir + '/sasa.xvg'
                elif os.path.exists(wdir + '/sasa.xvg'):
                    fname  = wdir + '/sasa.xvg'
                else:
                    print(fname, " does not exist! ")
                    continue

                # Open and parse file
                with open(fname) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines)

                yall = np.append(yall,data[:,1]) #append all y-data
                yavg += np.average(data[:,1])

            # append avg to yavg
            ysasa_avg = np.append(ysasa_avg,yavg/len(run_arr))
            # plot histogram across all cases
            sns.kdeplot(data=np.array(yall),label=solv_type,ax=ax1)

        # Save histogram plot
        plt.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_SASAdist.png',dpi=fig1.dpi)

        # Plot average
        fig2,ax2 = plt.subplots()
        ax2.set_ylabel(r'SASA (nm$^2$)')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(ysasa_avg),ax=ax2)
        change_width(ax2,0.5)
        fig2.savefig(fig_dir + '/'+biomass+'_SASAavg.png',dpi=fig2.dpi)

        if sasa == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                sasa_tot  = ysasa_avg[indx]
                fc_ss.write('%s\t%s\t%g\n' %(biomass,solv_type,\
                                             sasa_tot))
    if sasa == 2:
        fc_ss.close()
        df=pd.read_table(sasa_fyl)
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="SASA",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'SASA (nm$^2$)')
        change_width(axa,0.5)
        figa.savefig(fig_dir + '/'+'AllSASA.png',dpi=figa.dpi)
#------------------------------------------------------------------

# Plot hbonds
if hbonds != 0:

    print("Analyzing hydrogen bonds data")

    if hbonds == 2: # plot all averages together
        consol_dir = scr_dir + '/' + inp_type + '/consolidated'
        if not os.path.isdir(consol_dir):
            os.mkdir(consol_dir)

        # Write data and then read to concatenate
        hbc_fyl = consol_dir + '/hbdata.dat'
        fc_hb = open(hbc_fyl,'w')
        fc_hb.write('%s\t%s\t%s\t%s\t%s\n' %('Biomass','Solvent',\
                                             'HB_Wat','HB_Sol',\
                                             'HB_Tot'))
   
    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yhbw_avg = np.zeros(0); yhbs_avg = np.zeros(0)
        yhbt_avg = np.zeros(0)

        fig1, ax1 = plt.subplots()
        ax1.set_ylabel(r'$\langle$ #HB (Water) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        
        fig2, ax2 = plt.subplots()
        ax2.set_ylabel(r'$\langle$ #HB (Org. Solv.) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        
        fig3, ax3 = plt.subplots()
        ax3.set_ylabel(r'$\langle$ #HB (Total) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            y_wavg = 0; y_savg = 0; y_tavg = 0
            n_wfyl = 0; n_sfyl = 0; n_tfyl = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  solv_type)

                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file(s) for water hbond
                if glob.glob(anadir + '/wat_hbnum*') != []:
                    fwat_list = glob.glob(anadir + '/wat_hbnum*')
                elif glob.glob(wdir + '/wat_hbnum*') != []:
                    fwat_list = glob.glob(wdir + '/wat_hbnum*')
                else:
                    print("wat_hbnum does not exist!")

                for fname in fwat_list:
                    # Open and parse file
                    with open(fname) as fin:
                        lines = (line.lstrip() for line in fin \
                                 if not line.lstrip().startswith('#') and \
                                 not line.lstrip().startswith('@'))
                        data  = np.loadtxt(lines)

                    #append all hb-wat
                    y_wavg += np.sum(data[:,1]); n_wfyl += data[:,1].size
                    y_tavg += np.sum(data[:,1]); n_tfyl += data[:,1].size
                    
                # Check file(s) for solvent hbond
                if glob.glob(anadir + '/sol_hbnum*') != []:
                    fsol_list = glob.glob(anadir + '/sol_hbnum*')
                elif glob.glob(wdir + '/sol_hbnum*') != []:
                    fsol_list = glob.glob(wdir + '/sol_hbnum*')
                else:
                    print("sol_hbnum does not exist!")
                    continue


                for fname in fsol_list:
                    # Open and parse file
                    with open(fname) as fin:
                        lines = (line.lstrip() for line in fin \
                                 if not line.lstrip().startswith('#') and \
                                 not line.lstrip().startswith('@'))
                        data  = np.loadtxt(lines)

                    #append all hb-solv
                    y_savg += np.sum(data[:,1]); n_sfyl += data[:,1].size
                    y_tavg += np.sum(data[:,1])

                if n_wfyl != n_sfyl: # Sanity check
                    print(n_wfyl,n_sfyl)
                    print('total HB cannot be defined')
                    
            # append avg to overall averages
            yhbw_avg = np.append(yhbw_avg,y_wavg/n_wfyl)
            yhbs_avg = np.append(yhbs_avg,y_savg/n_sfyl)
            yhbt_avg = np.append(yhbt_avg,y_tavg/n_tfyl)

        # Plot averages
        sns.barplot(otyp_arr,np.array(yhbw_avg),ax=ax1)
        change_width(ax1,0.5)
        fig1.savefig(fig_dir + '/'+biomass+'_HBWat.png',dpi=fig1.dpi)

        sns.barplot(otyp_arr,np.array(yhbs_avg),ax=ax2)
        change_width(ax2,0.5)
        fig2.savefig(fig_dir + '/'+biomass+'_HBOrg.png',dpi=fig2.dpi)

        sns.barplot(otyp_arr,np.array(yhbt_avg),ax=ax3)
        change_width(ax3,0.5)
        fig3.savefig(fig_dir + '/'+biomass+'_HBTot.png',dpi=fig3.dpi)
        
        if hbonds == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                hbw = yhbw_avg[indx]; hbs = yhbs_avg[indx]
                hbt = yhbt_avg[indx]
                fc_hb.write('%s\t%s\t%g\t%g\t%g\n' %(biomass,solv_type,\
                                                     hbw,hbs,hbt))
    if hbonds == 2:
        fc_hb.close()
        df=pd.read_table(hbc_fyl)
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="HB_Tot",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'$\langle$ #HB (Total) $\rangle$')
        change_width(axa,0.5)
        figa.savefig(fig_dir + '/'+'AllHBTot.png',dpi=figa.dpi)
#------------------------------------------------------------------

# Plot RDF (bO4-solvent, bO4-water)
if rdf_bO4:
    
    print("Analyzing beta-O-4 RDF data")
    
    # Plot polymer-solvent RDF
    fig1,ax1 = plt.subplots()
    ax1.set_xlabel(r'$r$ (nm)')
    ax1.set_ylabel(r'$g_{\mathrm{pol-orgsol}}(r)$')
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    # Plot polymer-water RDF
    fig2,ax2 = plt.subplots()
    ax2.set_xlabel(r'$r$ (nm)')
    ax2.set_ylabel(r'$g_{\mathrm{pol-water}}(r)$')
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yrg_avg = np.zeros(0)

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            xdata = np.zeros(0)
            y1data = np.zeros(0); y2data = np.zeros(0)
            rdata = np.zeros(0) # for counting repeats

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  solv_type)
                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file
                if os.path.exists(anadir + '/rdfout2.xvg'):
                    fname  = anadir + '/rdfout2.xvg'
                elif os.path.exists(wdir + '/rdfout2.xvg'):
                    fname  = wdir + '/rdfout2.xvg'
                else:
                    print(fname, " does not exist! ")
                    continue

                # Open and parse file
                with open(fname) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines)

                if casenum == 0: #append all x/y-data
                    xdata = np.append(xdata,data[:,0])
                    y1data = np.append(y1data,data[:,1])
                    y2data = np.append(y2data,data[:,2])
                    omax_xval = np.amax(xdata)
                    olen_data = len(xdata)
                    rdata = np.full((len(xdata)),1,dtype=int)
                else:
                    if np.amax(data[:,0]) >= omax_xval:
                        for i_indx in range(olen_data):
                            y1data[i_indx] += data[i_indx,1]
                            y2data[i_indx] += data[i_indx,2]
                            rdata[i_indx] += 1
                        len_data = len(data[:,0])
                        max_xval = np.amax(data[:,0])
                        xdata = np.append(xdata,data[i_indx:len_data-1,0])
                        y1data = np.append(y1data,data[i_indx:len_data-1,1])
                        y2data = np.append(y2data,data[i_indx:len_data-1,1])
                        rnew  = np.full((1,len_data-olen_data),1,dtype=int)
                        rdata = np.append(rdata,rnew)
                        omax_xval = max_val
                        olen_data = len_data
                    else:
                        for i_indx in range(olen_data):
                            ydata[i_indx] += data[i_indx,1]
                            rdata[i_indx] += 1
                        
            # Divide ydata by rdata
            y1data /= rdata
            y2data /= rdata
            ax1.plot(xdata,y1data,label=solv_type)
            ax2.plot(xdata,y2data,label=solv_type)

        ax1.legend(loc=0)
        ax2.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_rdfsol.png',dpi=fig1.dpi)
        fig2.savefig(fig_dir + '/'+biomass+'_rdfwat.png',dpi=fig2.dpi)
#------------------------------------------------------------------

# Plot avg rg and rg distribution
if rad_gyr != 0:
    print("Analyzing Rg data")

    if rad_gyr == 2: # plot all averages together
        consol_dir = scr_dir + '/' + inp_type + '/consolidated'
        if not os.path.isdir(consol_dir):
            os.mkdir(consol_dir)

        # Write data and then read to concatenate
        rgall_fyl = consol_dir + '/Rgdata.dat'
        fc_rg = open(rgall_fyl,'w')
        fc_rg.write('%s\t%s\t%s\n' %('Biomass','Solvent','Rg'))

    fig1,ax1 = plt.subplots()
    ax1.set_xlabel(r'$r$' ' (nm)')
    ax1.set_ylabel(r'$R_g$' ' (nm)')
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()

    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yrg_avg = np.zeros(0)

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            yall = np.zeros(0)
            yavg = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  solv_type)
                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file
                if os.path.exists(anadir + '/rg_chain.xvg'):
                    fname  = anadir + '/rg_chain.xvg'
                elif os.path.exists(wdir + '/rg_chain.xvg'):
                    fname  = wdir + '/rg_chain.xvg'
                else:
                    print(fname, " does not exist! ")
                    continue

                # Open and parse file
                with open(fname) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines)

                yavg += np.average(data[:,1])
                yall = np.append(yall,data[:,1]) #append all y-data

            # append avg to yavg
            yrg_avg = np.append(yrg_avg,yavg/len(run_arr))

            # plot histogram across all cases
            sns.kdeplot(data=np.array(yall),label=solv_type,ax=ax1)

        # Save histogram plot
        plt.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_Rgdist.png',dpi=fig1.dpi)

        # Plot average
        fig2,ax2 = plt.subplots()
        ax2.set_ylabel(r'$R_g$ (nm)')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yrg_avg),ax=ax2)
        change_width(ax2,0.5)
        fig2.savefig(fig_dir + '/'+biomass+'_Rgavg.png',dpi=fig2.dpi)
        if sasa == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                rg_tot  = yrg_avg[indx]
                fc_rg.write('%s\t%s\t%g\n' %(biomass,solv_type,\
                                             rg_tot))
    if rad_gyr == 2:
        fc_rg.close()
        df=pd.read_table(rgall_fyl)
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="Rg",hue="Solvent",data=df,\
                    ax=axa)
        ax1.set_ylabel(r'$R_g$' ' (nm)')
        change_width(axa,0.5)
        figa.savefig(fig_dir + '/'+'AllRg.png',dpi=figa.dpi)
#------------------------------------------------------------------





