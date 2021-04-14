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
plt.rcParams.update({'font.size': 14}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
#plt.rcParams.update({'font.family': 'Times New Roman'})
#------------------------------------------------------------------

# Input data
inp_type  = 'cosolvents' # melts, solvents, cosolvents
disperse  = 'mono' # mono/poly; only for melts
biom_arr  = ['WT','COMT','MYB']#,'COMT','MYB'] # biomass type arr
otyp_arr  = ['EOH','GVL','THF']  # solvent arr for solvents/cosolvents
otyp_leg  = ['EtOH','GVL','THF']  # legend for solvents/cosolvents
run_arr   = [4,5,6] # number of independent runs for a given biomass
#------------------------------------------------------------------

# Plot keys 0,1,2 (0-None,1-separate,2-together)
sasa    = 0 # plot SASA distribution and average
rad_gyr = 0 # avg rg and rg distribution
rdf_pol = 1 # plot polymer-solvent/water RDF 0 or 1
rdf_bO4 = 0 # plot bO4-solvent/water RDF
hbonds  = 0 # hbonds of ALL-solv/water/total hbonds
hb_bo4  = 0 # hbonds of bo4-solv/water/total hbonds
hb_syr  = 0 # hbonds of syr-solv/water/total hbonds
hb_gua  = 0 # hbonds of gua-solv/water/total hbonds
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


    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        ysasa_avg = np.zeros(0)

        # Define axes labels for distribution
        fig1,ax1 = plt.subplots()
        ax1.set_xlabel(r'SASA (nm$^2$)')
        ax1.set_ylabel('Probability')
        plt.style.use('seaborn-colorblind')

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            yall = np.zeros(0)
            yavg = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'anafiles',solv_type)
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
            sns.kdeplot(data=np.array(yall),label=solv_leg,ax=ax1)

        # Save histogram plot
        ax1.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_SASAdist.png',dpi=fig1.dpi)
        plt.close(fig1)

        # Plot average
        fig2,ax2 = plt.subplots()
        ax2.set_ylabel(r'SASA (nm$^2$)')
        plt.style.use('seaborn-colorblind')
        sns.barplot(otyp_arr,np.array(ysasa_avg),ax=ax2)
        plt.tight_layout()
        change_width(ax2,0.5)
        fig2.savefig(fig_dir + '/'+biomass+'_SASAavg.png',dpi=fig2.dpi)
        plt.close(fig2)

        if sasa == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                solv_leg  = otyp_leg[indx]
                sasa_tot  = ysasa_avg[indx]
                fc_ss.write('%s\t%s\t%g\n' %(biomass,solv_leg,\
                                             sasa_tot))
    if sasa == 2:
        fc_ss.close()
        df=pd.read_table(sasa_fyl)
        maxval = df['SASA'].max()
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="SASA",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'SASA (nm$^2$)')
        axa.set_ylim([0,1.2*maxval])
        axa.legend(loc=0,ncol = len(axa.lines))
        change_width(axa,0.2)
        figa.savefig(fig_dir + '/'+'AllSASA.png',dpi=figa.dpi)
        plt.close(figa)
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

    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yrg_avg = np.zeros(0)

        # Define axes labels for distribution
        fig1,ax1 = plt.subplots()
        ax1.set_xlabel(r'$R_{g}$' ' (nm)')
        ax1.set_ylabel('Probability')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            yall = np.zeros(0)
            yavg = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'anafiles',solv_type)
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
            sns.kdeplot(data=np.array(yall),label=solv_leg,ax=ax1)

        # Save histogram plot
        ax1.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_Rgdist.png',dpi=fig1.dpi)
        plt.close(fig1)

        # Plot average
        fig2,ax2 = plt.subplots()
        ax2.set_ylabel(r'$R_g$ (nm)')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yrg_avg),ax=ax2)
        change_width(ax2,0.5)
        fig2.savefig(fig_dir + '/'+biomass+'_Rgavg.png',dpi=fig2.dpi)
        plt.close(fig2)

        if rad_gyr == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                solv_leg  = otyp_leg[indx]
                rg_tot  = yrg_avg[indx]
                fc_rg.write('%s\t%s\t%g\n' %(biomass,solv_leg,\
                                             rg_tot))
    if rad_gyr == 2:
        fc_rg.close()
        df=pd.read_table(rgall_fyl)
        maxval = df['Rg'].max()
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="Rg",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'$R_g$' ' (nm)')
        axa.set_ylim([0,1.2*maxval])
        axa.legend(loc=0,ncol = len(axa.lines))
        change_width(axa,0.2)
        figa.savefig(fig_dir + '/'+'AllRg.png',dpi=figa.dpi)
        plt.close(figa)
#------------------------------------------------------------------

# Plot RDF (bO4-solvent, bO4-water)
if rdf_pol:
    
    print("Analyzing whole polymer RDF data")
    
    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yrg_avg = np.zeros(0)

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

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            xdata = np.zeros(0)
            y1data = np.zeros(0); y2data = np.zeros(0)
            rdata = np.zeros(0) # for counting repeats

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'anafiles',solv_type)
                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file
                if os.path.exists(anadir + '/rdfout1.xvg'):
                    fname  = anadir + '/rdfout1.xvg'
                elif os.path.exists(wdir + '/rdfout1.xvg'):
                    fname  = wdir + '/rdfout1.xvg'
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
                        omax_xval = max_xval
                        olen_data = len_data
                    else:
                        for i_indx in range(len(data[:,0])):
                            y1data[i_indx] += data[i_indx,1]
                            y2data[i_indx] += data[i_indx,2]
                            rdata[i_indx] += 1
                        
            # Divide ydata by rdata
            y1data /= rdata
            y2data /= rdata
            ax1.plot(xdata,y1data,label=solv_leg)
            ax2.plot(xdata,y2data,label=solv_leg)

        ax1.legend(loc=0)
        ax2.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_rdfpsol.png',dpi=fig1.dpi)
        fig2.savefig(fig_dir + '/'+biomass+'_rdfpwat.png',dpi=fig2.dpi)
        plt.close(fig1)
        plt.close(fig2)
#------------------------------------------------------------------

# Plot RDF (bO4-solvent, bO4-water)
if rdf_bO4:
    
    print("Analyzing beta-O-4 RDF data")
    
    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yrg_avg = np.zeros(0)

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

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            xdata = np.zeros(0)
            y1data = np.zeros(0); y2data = np.zeros(0)
            rdata = np.zeros(0) # for counting repeats

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'anafiles',solv_type)
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
                        omax_xval = max_xval
                        olen_data = len_data
                    else:
                        for i_indx in range(len(data[:,0])):
                            y1data[i_indx] += data[i_indx,1]
                            y2data[i_indx] += data[i_indx,2]
                            rdata[i_indx] += 1
                        
            # Divide ydata by rdata
            y1data /= rdata
            y2data /= rdata
            ax1.plot(xdata,y1data,label=solv_leg)
            ax2.plot(xdata,y2data,label=solv_leg)

        ax1.legend(loc=0)
        ax2.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_rdfsol.png',dpi=fig1.dpi)
        fig2.savefig(fig_dir + '/'+biomass+'_rdfwat.png',dpi=fig2.dpi)
        plt.close(fig1)
        plt.close(fig2)
#------------------------------------------------------------------

# Plot hbonds
if hbonds != 0:

    print("Analyzing total hydrogen bonds data")

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

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            y_wavg = 0; y_savg = 0; y_tavg = 0
            n_wfyl = 0; n_sfyl = 0; n_tfyl = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'anafiles',solv_type)

                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file(s) for water hbond
                if glob.glob(anadir + '/wat_hbnum*') != []:
                    fwat_list = glob.glob(anadir + '/wat_hbnum*')
                elif glob.glob(wdir + '/wat_hbnum*') != []:
                    fwat_list = glob.glob(wdir + '/wat_hbnum*')
                else:
                    print("wat_hbnum does not exist!")
                    continue

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
        fig1, ax1 = plt.subplots()
        ax1.set_ylabel(r'$\langle$ #HB (Water) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yhbw_avg),ax=ax1)
        change_width(ax1,0.5)
        fig1.savefig(fig_dir + '/'+biomass+'_HBWat.png',dpi=fig1.dpi)
        
        fig2, ax2 = plt.subplots()
        ax2.set_ylabel(r'#HBs with the organic solvent')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yhbs_avg),ax=ax2)
        change_width(ax2,0.5)
        fig2.savefig(fig_dir + '/'+biomass+'_HBOrg.png',dpi=fig2.dpi)
        
        fig3, ax3 = plt.subplots()
        ax3.set_ylabel(r'$\langle$ #HB (Total) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yhbt_avg),ax=ax3)
        change_width(ax3,0.5)
        fig3.savefig(fig_dir + '/'+biomass+'_HBTot.png',dpi=fig3.dpi)

        plt.close(fig1); plt.close(fig2); plt.close(fig3)
        # Save average data        
        if hbonds == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                solv_leg  = otyp_leg[indx]
                hbw = yhbw_avg[indx]; hbs = yhbs_avg[indx]
                hbt = yhbt_avg[indx]
                fc_hb.write('%s\t%s\t%g\t%g\t%g\n' %(biomass,solv_leg,\
                                                     hbw,hbs,hbt))
    if hbonds == 2:
        fc_hb.close()
        df=pd.read_table(hbc_fyl)
        maxval = df['HB_Tot'].max()
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="HB_Tot",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'$\langle$ #HB (Total) $\rangle$')
        axa.set_ylim([0,1.2*maxval])
        axa.legend(loc=0,ncol = len(axa.lines))
        change_width(axa,0.2)
        figa.savefig(fig_dir + '/'+'AllHBTot.png',dpi=figa.dpi)

        maxval = df['HB_Wat'].max()
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="HB_Wat",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'$\langle$ #Water-HB $\rangle$')
        axa.set_ylim([0,1.2*maxval])
        axa.legend(loc=0,ncol = len(axa.lines))
        change_width(axa,0.2)
        figa.savefig(fig_dir + '/'+'AllHBWat.png',dpi=figa.dpi)

        maxval = df['HB_Sol'].max()
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="HB_Sol",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'$\langle$ #HBs with the Organic Solvent $\rangle$')
        axa.set_ylim([0,1.2*maxval])
        axa.set_xlabel('')
        plt.minorticks_on()
        axa.tick_params(axis='x', which='minor', bottom=False)
        axa.legend(loc=0,ncol = len(axa.lines))
        change_width(axa,0.2)
        figa.savefig(fig_dir + '/'+'AllHBOsol.png',dpi=figa.dpi)
#------------------------------------------------------------------

# Plot beta-O-4 hydrogen bonds with water and solvent
if hb_bo4 != 0:

    print("Analyzing beta-O-4 hydrogen bonds data")

    if hb_bo4 == 2: # plot all averages together
        consol_dir = scr_dir + '/' + inp_type + '/consolidated'
        if not os.path.isdir(consol_dir):
            os.mkdir(consol_dir)

        # Write data and then read to concatenate
        hbc_fyl = consol_dir + '/hb_bo4.dat'
        fc_hb = open(hbc_fyl,'w')
        fc_hb.write('%s\t%s\t%s\t%s\t%s\n' %('Biomass','Solvent',\
                                             'HB_BWat','HB_BSol',\
                                             'HB_BTot'))
   
    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yhbw_avg = np.zeros(0); yhbs_avg = np.zeros(0)
        yhbt_avg = np.zeros(0)

        for sol_indx in range(len(otyp_arr)): # loop in solvents
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            y_wavg = 0; y_savg = 0; y_tavg = 0
            n_wfyl = 0; n_sfyl = 0; n_tfyl = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'hbfiles',solv_type)

                print("Analyzing: ", biomass,solv_type,run_arr[casenum])

                # Check file(s) for water hbond
                if glob.glob(anadir + '/wat_hbbo4?.xvg') != []:
                    fwat_list = glob.glob(anadir + '/wat_hbbo4?.xvg')
                elif glob.glob(wdir + '/wat_hbbo4?.xvg') != []:
                    fwat_list = glob.glob(wdir + '/wat_hbbo4?.xvg')
                else:
                    print("wat_hbbo4 does not exist!")
                    continue

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
                if glob.glob(anadir + '/sol_hbbo4?.xvg') != []:
                    fsol_list = glob.glob(anadir + '/sol_hbbo4?.xvg')
                elif glob.glob(wdir + '/sol_hbbo4?.xvg') != []:
                    fsol_list = glob.glob(wdir + '/sol_hbbo4?.xvg')
                else:
                    print("sol_hbbo4 does not exist!")
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

            # append avg to overall averages
            if n_sfyl == 0 and y_savg == 0:
                print("WARNING: NO SOLVENT FILES WERE FOUND")
                yhbs_avg = np.append(yhbs_avg,0)
            elif n_wfyl != n_sfyl: # Sanity check
                print(n_wfyl,n_sfyl)
                print('total HB for beta-O-4 cannot be defined')
            else:
                yhbs_avg = np.append(yhbs_avg,y_savg/n_sfyl)

            yhbt_avg = np.append(yhbt_avg,y_tavg/n_tfyl)
            yhbw_avg = np.append(yhbw_avg,y_wavg/n_wfyl)

        # Plot averages
        fig1, ax1 = plt.subplots()
        ax1.set_ylabel(r'$\langle$ #HB ($\beta$-O-4 - Water) $\rangle$')
        plt.style.use('seaborn-colorblind')
        sns.barplot(otyp_arr,np.array(yhbw_avg),ax=ax1)
        plt.tight_layout()
        change_width(ax1,0.5)
        fig1.savefig(fig_dir + '/'+biomass+'_HBbo4Wat.png',dpi=fig1.dpi)
               
        if not all(ysv == 0 for ysv in yhbs_avg):
            fig2, ax2 = plt.subplots()
            ax2.set_ylabel(r'$\langle$ #HB ($\beta$-O-4 - Org. Solv.) $\rangle$')
            plt.style.use('seaborn-colorblind')
            sns.barplot(otyp_arr,np.array(yhbs_avg),ax=ax2)
            plt.tight_layout()
            change_width(ax2,0.5)
            fig2.savefig(fig_dir + '/'+biomass+'_HBbo4Org.png',dpi=fig2.dpi)
            plt.close(fig2)

        fig3, ax3 = plt.subplots()
        ax3.set_ylabel(r'$\langle$ #HB ($\beta$-O-4 - Total) $\rangle$')
        plt.style.use('seaborn-colorblind')
        sns.barplot(otyp_arr,np.array(yhbt_avg),ax=ax3)
        plt.tight_layout()
        change_width(ax3,0.5)
        fig3.savefig(fig_dir + '/'+biomass+'_HBbo4Tot.png',dpi=fig3.dpi)

        plt.close(fig1); plt.close(fig3)        
        # Save average data
        if hb_bo4 == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                solv_leg  = otyp_leg[indx]
                hbw = yhbw_avg[indx]; hbs = yhbs_avg[indx]
                hbt = yhbt_avg[indx]
                fc_hb.write('%s\t%s\t%g\t%g\t%g\n' %(biomass,solv_leg,\
                                                     hbw,hbs,hbt))
    if hb_bo4 == 2:
        fc_hb.close()
        df=pd.read_table(hbc_fyl)
        maxval = df['HB_BTot'].max()
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        sns.barplot(x="Biomass",y="HB_BTot",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'$\langle$ #HB ($\beta$-O-4) $\rangle$')
        axa.set_ylim([0,1.2*maxval])
        axa.legend(loc=0,ncol = len(axa.lines))
        plt.tight_layout()
        change_width(axa,0.2)
        figa.savefig(fig_dir + '/'+'AllHB_bo4Tot.png',dpi=figa.dpi)
        plt.close(figa)
#------------------------------------------------------------------

# Plot syringyl hydrogen bonds with water and solvent
if hb_syr != 0:

    print("Analyzing hydrogen bonds with SYR monomers")

    if hb_syr == 2: # plot all averages together
        consol_dir = scr_dir + '/' + inp_type + '/consolidated'
        if not os.path.isdir(consol_dir):
            os.mkdir(consol_dir)

        # Write data and then read to concatenate
        hbc_fyl = consol_dir + '/hb_syr.dat'
        fc_hb = open(hbc_fyl,'w')
        fc_hb.write('%s\t%s\t%s\t%s\t%s\n' %('Biomass','Solvent',\
                                             'HB_SWat','HB_SSol',\
                                             'HB_STot'))
   
    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yhbw_avg = np.zeros(0); yhbs_avg = np.zeros(0)
        yhbt_avg = np.zeros(0)

        for sol_indx in range(len(otyp_arr)): # loop in solvents (SYR)
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            y_wavg = 0; y_savg = 0; y_tavg = 0
            n_wfyl = 0; n_sfyl = 0; n_tfyl = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'hbfiles',solv_type)

                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file(s) for SYR-water hbond
                if glob.glob(anadir + '/wat_syrmon?.xvg') != []:
                    fwat_list = glob.glob(anadir + '/wat_syrmon?.xvg')
                elif glob.glob(wdir + '/wat_syrmon?.xvg') != []:
                    fwat_list = glob.glob(wdir + '/wat_syrmon?.xvg')
                else:
                    print("wat_syrmon does not exist!")
                    continue

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
                if glob.glob(anadir + '/sol_syrmon?.xvg') != []:
                    fsol_list = glob.glob(anadir + '/sol_syrmon?.xvg')
                elif glob.glob(wdir + '/sol_syrmon?.xvg') != []:
                    fsol_list = glob.glob(wdir + '/sol_syrmon?.xvg')
                else:
                    print("sol_syrmon does not exist!")
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
                
                print(casenum, n_wfyl, n_sfyl)
            # append avg to overall averages
            if n_sfyl == 0 and y_savg == 0:
                print("WARNING: NO SOLVENT FILES WERE FOUND")
                yhbs_avg = np.append(yhbs_avg,0)
            elif n_wfyl != n_sfyl: # Sanity check
                print(n_wfyl,n_sfyl)
                print('total HB for SYR cannot be defined')
            else:
                yhbs_avg = np.append(yhbs_avg,y_savg/n_sfyl)
                
            yhbt_avg = np.append(yhbt_avg,y_tavg/n_tfyl)
            yhbw_avg = np.append(yhbw_avg,y_wavg/n_wfyl)


        # Plot averages
        fig1, ax1 = plt.subplots()
        ax1.set_ylabel(r'$\langle$ #HB (SYR-Water) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yhbw_avg),ax=ax1)
        change_width(ax1,0.5)
        fig1.savefig(fig_dir + '/'+biomass+'_HBSyrWat.png',dpi=fig1.dpi)

        if not all(ysv == 0 for ysv in yhbs_avg):        
            fig2, ax2 = plt.subplots()
            ax2.set_ylabel(r'$\langle$ #HB (SYR-Org. Solv.) $\rangle$')
            plt.style.use('seaborn-colorblind')
            plt.tight_layout()
            sns.barplot(otyp_arr,np.array(yhbs_avg),ax=ax2)
            change_width(ax2,0.5)
            fig2.savefig(fig_dir + '/'+biomass+'_HBSyrOrg.png',dpi=fig2.dpi)
            plt.close(fig2)        

        fig3, ax3 = plt.subplots()
        ax3.set_ylabel(r'$\langle$ #HB (SYR-Total) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yhbt_avg),ax=ax3)
        change_width(ax3,0.5)
        fig3.savefig(fig_dir + '/'+biomass+'_HBSyrTot.png',dpi=fig3.dpi)
        plt.close(fig1); plt.close(fig3)        

        # Save average data        
        if hb_syr == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                solv_leg  = otyp_leg[indx]
                hbw = yhbw_avg[indx]; hbs = yhbs_avg[indx]
                hbt = yhbt_avg[indx]
                fc_hb.write('%s\t%s\t%g\t%g\t%g\n' %(biomass,solv_leg,\
                                                     hbw,hbs,hbt))
    if hb_syr == 2:
        fc_hb.close()
        df=pd.read_table(hbc_fyl)
        maxval = df['HB_STot'].max()
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="HB_STot",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'$\langle$ #HB (SYR) $\rangle$')
        axa.set_ylim([0,1.2*maxval])
        axa.legend(loc=0,ncol = len(axa.lines))
        change_width(axa,0.2)
        figa.savefig(fig_dir + '/'+'AllHB_SyrTot.png',dpi=figa.dpi)
        plt.close(figa)
#------------------------------------------------------------------

# Plot guaicyl hydrogen bonds with water and solvent
if hb_gua != 0:

    print("Analyzing hydrogen bonds with GUA monomers")

    if hb_gua == 2: # plot all averages together
        consol_dir = scr_dir + '/' + inp_type + '/consolidated'
        if not os.path.isdir(consol_dir):
            os.mkdir(consol_dir)

        # Write data and then read to concatenate
        hbc_fyl = consol_dir + '/hb_gua.dat'
        fc_hb = open(hbc_fyl,'w')
        fc_hb.write('%s\t%s\t%s\t%s\t%s\n' %('Biomass','Solvent',\
                                             'HB_GWat','HB_GSol',\
                                             'HB_GTot'))
   
    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yhbw_avg = np.zeros(0); yhbs_avg = np.zeros(0)
        yhbt_avg = np.zeros(0)

        for sol_indx in range(len(otyp_arr)): # loop in solvents (SYR)
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            y_wavg = 0; y_savg = 0; y_tavg = 0
            n_wfyl = 0; n_sfyl = 0; n_tfyl = 0

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'hbfiles',solv_type)

                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file(s) for GUA-water hbond
                if glob.glob(anadir + '/wat_guamon?.xvg') != []:
                    fwat_list = glob.glob(anadir + '/wat_guamon?.xvg')
                elif glob.glob(wdir + '/wat_guamon?.xvg') != []:
                    fwat_list = glob.glob(wdir + '/wat_guamon?.xvg')
                else:
                    print("wat_guamon does not exist!")
                    exist

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
                if glob.glob(anadir + '/sol_guamon?.xvg') != []:
                    fsol_list = glob.glob(anadir + '/sol_syrmon?.xvg')
                elif glob.glob(wdir + '/sol_guamon?.xvg') != []:
                    fsol_list = glob.glob(wdir + '/sol_guamon?.xvg')
                else:
                    print("sol_guamon does not exist!")
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

                    
            # append avg to overall averages
            if n_sfyl == 0 and y_savg == 0:
                print("WARNING: NO SOLVENT FILES WERE FOUND")
                yhbs_avg = np.append(yhbs_avg,0)
            elif n_wfyl != n_sfyl: # Sanity check
                print(n_wfyl,n_sfyl)
                print('total HB for GUA cannot be defined')
            else:
                yhbs_avg = np.append(yhbs_avg,y_savg/n_sfyl)
                
            yhbt_avg = np.append(yhbt_avg,y_tavg/n_tfyl)
            yhbw_avg = np.append(yhbw_avg,y_wavg/n_wfyl)

        # Plot averages
        fig1, ax1 = plt.subplots()
        ax1.set_ylabel(r'$\langle$ #HB (GUA-Water) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yhbw_avg),ax=ax1)
        change_width(ax1,0.5)
        fig1.savefig(fig_dir + '/'+biomass+'_HBGuaWat.png',dpi=fig1.dpi)
        
        if not all(ysv == 0 for ysv in yhbs_avg):        
            fig2, ax2 = plt.subplots()
            ax2.set_ylabel(r'$\langle$ #HB (GUA-Org. Solv.) $\rangle$')
            plt.style.use('seaborn-colorblind')
            plt.tight_layout()
            sns.barplot(otyp_arr,np.array(yhbs_avg),ax=ax2)
            change_width(ax2,0.5)
            fig2.savefig(fig_dir + '/'+biomass+'_HBGuaOrg.png',dpi=fig2.dpi)
            plt.close(fig2)

        fig3, ax3 = plt.subplots()
        ax3.set_ylabel(r'$\langle$ #HB (GUA-Total) $\rangle$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(otyp_arr,np.array(yhbt_avg),ax=ax3)
        change_width(ax3,0.5)
        fig3.savefig(fig_dir + '/'+biomass+'_HBGuaTot.png',dpi=fig3.dpi)
        plt.close(fig1); plt.close(fig3)        

        # Save average data        
        if hb_gua == 2:
            for indx in range(len(otyp_arr)): # loop in solvents
                solv_type = otyp_arr[indx]
                solv_leg  = otyp_leg[indx]
                hbw = yhbw_avg[indx]; hbs = yhbs_avg[indx]
                hbt = yhbt_avg[indx]
                fc_hb.write('%s\t%s\t%g\t%g\t%g\n' %(biomass,solv_leg,\
                                                     hbw,hbs,hbt))
    if hb_gua == 2:
        fc_hb.close()
        df=pd.read_table(hbc_fyl)
        maxval = df['HB_GTot'].max()
        figa, axa = plt.subplots()
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()
        sns.barplot(x="Biomass",y="HB_GTot",hue="Solvent",data=df,\
                    ax=axa)
        axa.set_ylabel(r'$\langle$ #HB (GUA) $\rangle$')
        axa.set_ylim([0,1.2*maxval])
        axa.legend(loc=0,ncol = len(axa.lines))
        change_width(axa,0.2)
        figa.savefig(fig_dir + '/'+'AllHB_GuaTot.png',dpi=figa.dpi)
        plt.close(figa)
#------------------------------------------------------------------
