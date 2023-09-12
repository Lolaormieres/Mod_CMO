#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:30:16 2023

@author: ormieresl
"""
#! /usr/bin/python

# coding: utf-8

#import bronx
from vortex import toolbox
import common
import olive
import common.util.usepygram
import usevortex
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import epygram
from datetime import datetime, timedelta
import vortex
import pandas as pd
import os
import numpy.ma as ma
epygram.init_env()
from matplotlib import font_manager as fm, rcParams




# name=['SFX.SIC','SFX.SST','SFX.ICEHSI_1','SFX.ICETSF_1']
# MINMAX=[0,1]
# MINMAX=[270,280] #rajout
# MINMAX=[[0,1],[270,280],[0,4.],[-40,1]]
# MINMAX=[[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]]

#name = ['SFX.SST']
#name_mer = ['SURFSST.CLIM.']
name = ['SFX.TEMP_OC1', 'SFX.TEMP_OC2', 'SFX.TEMP_OC3', 'SFX.TEMP_OC4', 'SFX.TEMP_OC5', 'SFX.TEMP_OC6', 'SFX.TEMP_OC7', 'SFX.TEMP_OC8', 'SFX.TEMP_OC9', 'SFX.TEMP_OC10', 'SFX.TEMP_OC11', 'SFX.TEMP_OC12']#,
        #'SFX.TEMP_OC13','SFX.TEMP_OC14', 'SFX.TEMP_OC15', 'SFX.TEMP_OC16', 'SFX.TEMP_OC17', 'SFX.TEMP_OC18', 'SFX.TEMP_OC19', 'SFX.TEMP_OC20', 'SFX.TEMP_OC21', 'SFX.TEMP_OC22', 'SFX.TEMP_OC23', 'SFX.TEMP_OC24', 'SFX.TEMP_OC25', 'SFX.TEMP_OC26']
name_sal = ['SFX.SALT_OC1', 'SFX.SALT_OC2', 'SFX.SALT_OC3', 'SFX.SALT_OC4', 'SFX.SALT_OC5', 'SFX.SALT_OC6', 'SFX.SALT_OC7', 'SFX.SALT_OC8', 'SFX.SALT_OC9', 'SFX.SALT_OC10', 'SFX.SALT_OC11', 'SFX.SALT_OC12', 'SFX.SALT_OC13',
            'SFX.SALT_OC14', 'SFX.SALT_OC15', 'SFX.SALT_OC16', 'SFX.SALT_OC17', 'SFX.SALT_OC18', 'SFX.SALT_OC19', 'SFX.SALT_OC20', 'SFX.SALT_OC21', 'SFX.SALT_OC22', 'SFX.SALT_OC23', 'SFX.SALT_OC24', 'SFX.SALT_OC25', 'SFX.SALT_OC26']
name_mer = ['SURF.THETAO1']#, 'SURF.THETAO2', 'SURF.THETAO3', 'SURF.THETAO4', 'SURF.THETAO5', 'SURF.THETAO6', 'SURF.THETAO7', 'SURF.THETAO8', 'SURF.THETAO9', 'SURF.THETAO10', 'SURF.THETAO11', 'SURF.THETAO12']#,
            # 'SURF.THETAO13','SURF.THETAO14', 'SURF.THETAO15', 'SURF.THETAO16', 'SURF.THETAO17', 'SURF.THETAO18', 'SURF.THETAO19', 'SURF.THETAO20', 'SURF.THETAO21', 'SURF.THETAO22', 'SURF.THETAO23', 'SURF.THETAO24', 'SURF.THETAO25', 'SURF.THETAO26']
name_sal_mer = ['SURF.SALINO1', 'SURF.SALINO2', 'SURF.SALINO3', 'SURF.SALINO4', 'SURF.SALINO5', 'SURF.SALINO6', 'SURF.SALINO7', 'SURF.SALINO8', 'SURF.SALINO9', 'SURF.SALINO10', 'SURF.SALINO11', 'SURF.SALINO12', 'SURF.SALINO13',
                'SURF.SALINO14', 'SURF.SALINO15', 'SURF.SALINO16', 'SURF.SALINO17', 'SURF.SALINO18', 'SURF.SALINO19', 'SURF.SALINO20', 'SURF.SALINO21', 'SURF.SALINO22', 'SURF.SALINO23', 'SURF.SALINO24', 'SURF.SALINO25', 'SURF.SALINO26']
name_u = ['SFX.UCUR_OC1', 'SFX.UCUR_OC2', 'SFX.UCUR_OC3', 'SFX.UCUR_OC4', 'SFX.UCUR_OC5', 'SFX.UCUR_OC6', 'SFX.UCUR_OC7', 'SFX.UCUR_OC8', 'SFX.UCUR_OC9', 'SFX.UCUR_OC10', 'SFX.UCUR_OC11', 'SFX.UCUR_OC12', 'SFX.UCUR_OC13',
          'SFX.UCUR_OC14', 'SFX.UCUR_OC15', 'SFX.UCUR_OC16', 'SFX.UCUR_OC17', 'SFX.UCUR_OC18', 'SFX.UCUR_OC19', 'SFX.UCUR_OC20', 'SFX.UCUR_OC21', 'SFX.UCUR_OC22', 'SFX.UCUR_OC23', 'SFX.UCUR_OC24', 'SFX.UCUR_OC25', 'SFX.UCUR_OC26']
name_v = ['SFX.VCUR_OC1', 'SFX.VCUR_OC2', 'SFX.VCUR_OC3', 'SFX.VCUR_OC4', 'SFX.VCUR_OC5', 'SFX.VCUR_OC6', 'SFX.VCUR_OC7', 'SFX.VCUR_OC8', 'SFX.VCUR_OC9', 'SFX.VCUR_OC10', 'SFX.VCUR_OC11', 'SFX.VCUR_OC12', 'SFX.VCUR_OC13',
          'SFX.VCUR_OC14', 'SFX.VCUR_OC15', 'SFX.VCUR_OC16', 'SFX.VCUR_OC17', 'SFX.VCUR_OC18', 'SFX.VCUR_OC19', 'SFX.VCUR_OC20', 'SFX.VCUR_OC21', 'SFX.VCUR_OC22', 'SFX.VCUR_OC23', 'SFX.VCUR_OC24', 'SFX.VCUR_OC25', 'SFX.VCUR_OC26']


## Lecture des niveaux de la CMO
filename = 'level2.txt'
level = np.loadtxt(filename, delimiter=',', dtype=float)
print('level', level)
print(type(level))
print(np.shape(level))
level_12 = level[0:12]
print(level_12)
level_26 = level[0:26]



## ARP - assimilation
cfields = ''
cblock = 'surfan'
ckind = 'analysis'
cut = 'assim'
# cblock='forecast'
# ckind='historic'
## ARP - production
cblock_for = 'forecast'
ckind_for = 'historic'
cmodel = 'surfex'
cterm = '0'
cmodel_for = 'surfex'
cut_for = 'production'
cterm_for = '96'
ech=''
## Mercator 
cblock_oc = 'c933'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
#cfields_oc = 'sst'
cmodel_oc = 'arpege'
cut_oc = 'assim'

# geometry='global798c22'
geometry = 'global1798'
# chemin='/d0/images/ormieresl/'
ech=cterm

dir_carole='/home/ormieresl/Routines/gmapfactor/'
dir_lola= '/home/ormieresl/Routines/DF_TEMP_mapfactor_mask/'

dir_fig1='/cnrm/recyf/Data/users/ormieresl/plot_score_water_column_temp/'
dir_fig2='/cnrm/recyf/Data/users/ormieresl/plot_score_water_column_temp_1month/'
dir_fig2='/cnrm/recyf/Data/users/ormieresl/plot_score_water_column_temp_septembre_vs_fevrier/'

## Definition color
color_GL55='#648FFF'  ## GL55
color_GKOR='#785EF0'  ## GKOR
color_GKPH='red'  ## GKPH '#DC267F'
color_GM6E='#FE6100'  ## GM6E
color_GMOT='#009E73'  ## GMOT
color_GK7C='dimgrey' 
color_GN3C='gold'
color_GN84='lightseagreen'
color_GNBA='midnightblue'

#%%
# =============================================================================
# ## Lectur ficher fa
# =============================================================================

def get_file(xpid, date, rep, typfic, fields, mod, chain, ech):
    # MODEL
    fp = {'experiment': xpid, 'kind': typfic, 'block': rep, 'fields': fields,
          'date': date, 'term': ech, 'geometry': geometry, 'format': 'fa',
          'local': 'tmpech.fic', 'cutoff': chain, 'vapp': 'arpege',
          'vconf': '4dvarfr', 'model': mod, 'namespace': 'vortex.multi.fr', 'filling': 'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpech.fic', openmode='r')
    os.remove('tmpech.fic')
    #return fcref.contents.data
    #fcref.container.filename
    return r

    
#%%
zone=['eurat','med','tropiques','nordat','hn20','hs20','glob']
zone_name=['Eurat','Meditteranean','Tropics', 'North Atlantic', 'North Hemisphere', 'South Hemisphere', 'Global' ]
#Period Temperatures
#bar_tick=[np.linspace(17, 18, 11),np.linspace(23.7, 24.9,7),np.linspace(26.6, 27.1, 11), np.linspace(18.6, 19.6, 11), np.linspace(16, 17, 11), np.linspace(10.5, 11.5, 11),np.linspace(18, 19, 11)]
# bar_tick=[np.linspace(17.6, 18, 9),np.linspace(23, 23.6,13),np.linspace(26.8, 27.2, 9),np.linspace(18.4, 18.9, 6),np.linspace(15.8, 16.3, 6),np.linspace(11, 11.5, 6),np.linspace(18.4, 18.7, 7)] # Octobre
# bar_tick=[np.linspace(19, 19.7, 8),np.linspace(24.5, 26.2, 11), np.linspace(26.5, 27, 6),np.linspace(19.5, 20.3, 9),np.linspace(17.1, 17.7, 7),np.linspace(10.8, 11.1, 7),np.linspace(18.7, 19, 7)] # septembre
# bar_tick=[np.linspace(12.4, 12.9, 11),np.linspace(14.7, 15.1, 11), np.linspace(27, 27.5, 6),np.linspace(13.7, 14.2, 11),np.linspace(10.9, 11.3, 9),np.linspace(14.7, 15.2, 11),np.linspace(17.7, 18.1, 9)] # fevrier

#Domain temperature
#,'Tropics','South Hem.','Global','Med']
# bar_tick=np.linspace(23, 25, 11) # Med
# bar_ticks=np.linspace(27, 28, 11) # Tropics
# bar_tick=np.linspace(18, 19, 11) # Eurat

expe=['GL55','GKOR','GKPH','GN84']
expe=['GN84']

mois='fevrier'
period='fevrier_mars'

mois='septembre'
period='septembre_octobre'

mois='1day_02092022'
period='septembre_octobre'
#%% TEMP EXPE
for izone in range(len(zone)):
    print(izone, zone[izone])
    
    ## FIGURE
    plt.rcParams["figure.figsize"] = [5, 4.50]
    plt.rcParams["figure.autolayout"] = True
    fig, axs = plt.subplots(nrows = 1,
                                    ncols = 1, 
                                    sharex = False, 
                                    figsize=(15, 7),dpi=150)


    df_ech = pd.read_csv(dir_lola+str(expe[0])+'_'+str(zone[izone])+'_ech_'+str(period)+'_prev-mer'+'.csv')
    df_ech.head()
    df_ech.dtypes
    #df_ech.columns.values GET name of columns
    df_ech=df_ech.loc[((df_ech['expe']==expe[0]) )][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar','temp_exp_mean','temp_merc_mean']]
    # df_ech=df_ech[5915:12239] # octobre
    df_ech=df_ech[0:5915] # septembre
    df_ech=df_ech[0:204] # 1day
    
    # DROP first raw
    # df_ech.drop(index=df_ech.index[0], axis=0, inplace=True)
    for niv in name:
        # print(niv)
        ech2 = np.unique(df_ech['ech'])  ## fct get only one time repeat variables
        N_level = name
        
        ## CREATE DF fill up with naans 3 cols (17 raws)
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'temp_exp_mean': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res_cmo = pd.DataFrame(data)
        
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'Tmer': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res_mer = pd.DataFrame(data)
    
        # MEAN each due date
        i = 0
        for e in ech2:
            for N in N_level:
                cmo_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['var']==N)].mean()['temp_exp_mean'])
                # print(biais_mean)
                df_res_cmo.iloc[i] = np.array([N,e,round(cmo_mean,5)])
                i = i + 1
                
                
        # MEAN each due date
        i = 0
        for e in ech2:
            for N in N_level:
                temp_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['var']==N)].mean()['temp_merc_mean'])
                # print(biais_mean)
                df_res_mer.iloc[i] = np.array([N,e,round(temp_mean,5)])
                i = i + 1

        # FIG
        df_res_cmo["temp_exp_mean"]=df_res_cmo["temp_exp_mean"].astype(float)
       
        N=1
        ## PLOT CONTOUR
        ax=axs
        ## DATE - XTICK
        ech=len(ech2)
        L=len(ech2)
        lst = list(np.arange(1,L+1))
        xx=lst
        list_ech=np.arange(0,102,6)
        xind = xx[0:len(xx):N] # Pas des dates plot
        xlabels = list_ech[0:len(list_ech):N]
        x = lst
        y = level[0:12]
        
        X, Y = np.meshgrid(x, y)
        
        temp=[]
        temp=df_res_cmo.temp_exp_mean
        temp=np.array(temp)
        
        ## Reshape bias en 2D
        reshaped_temp = temp.reshape(ech,12).T
        Z = reshaped_temp
        
        ## PLOT PROPRE
        ax.set_title('Zone: '+str(zone_name[izone]),fontsize=30,y=1.06)
        ax.set_xlabel('Hours',fontsize=30)
        ax.set_xticks(lst,
                labels=xlabels)
        axs.set_ylabel('Depth (m)',fontsize=30)
    
        # cbar_ax = fig.add_axes([0.3, -0.17, 0.5, 0.08])  # [left,bottom,width,height]
        cbar_ax = fig.add_axes([1.05, 0.095, 0.02, 0.78])
      
        #CS1=ax.contourf(X, Y, Z,cmap='rainbow', levels=bar_tick[izone], linewidths=1)
        
        CS1=ax.contourf(X, Y, Z,cmap='rainbow', linewidths=1)
        clb=fig.colorbar(CS1, spacing='uniform', cax=cbar_ax, orientation='vertical')
        cbar_ax.tick_params(labelsize=26)
        ax.tick_params(axis='y', labelsize=26)
        ax.tick_params(axis='x', labelsize=26)
        clb.ax.set_title('Temperature (°C)',fontsize=28, y=1.05)
       
        # plt.legend()
        fig.colorbar(CS1, cax=cbar_ax, orientation='vertical')
        # fig.suptitle('Evolution du biais de température par rapport à la Reference Mercator -'+str(expe), horizontalalignment = 'center',fontsize=15)
        fig.suptitle(str(mois)+'- averaged temperature in 1d model - '+str(expe), horizontalalignment = 'center',fontsize=28, y=1.02)
    print('min',min(temp))
    print('max',max(temp))
    print('averaged',np.mean(temp))
    
    figname = dir_fig2+str(expe[0])+'_xptemp_column_ech_HN_'+str(zone[izone])+'_'+str(mois)+'.png'
    fig.savefig(figname,dpi=150, format='png',bbox_inches='tight')
    figname = dir_fig2+str(expe[0])+'_xptemp_column_ech_HN_'+str(zone[izone])+'_'+str(mois)+'.pdf'
    fig.savefig(figname,dpi=150, format='pdf',bbox_inches='tight')

#%%
# zone=['eurat']#,'med','tropiques','nordat','hn20','hs20','glob']
# zone_name=['Eurat']#,'Meditteranean','Tropics', 'North Atlantic', 'North Hemisphere', 'South Hemisphere', 'Global' ]
#TEMP MERCATOR
# bar_tick=[np.linspace(17.6, 18, 9)]
#bar_tick=[np.linspace(18, 19, 11),np.linspace(23.7, 24.9,7),np.linspace(26.6, 27.1, 11), np.linspace(18.6,19.6, 11), np.linspace(16, 17, 11), np.linspace(10.5, 11.5, 11),np.linspace(18, 19, 11)]
#bar_tick=[np.linspace(17.6, 18, 9),np.linspace(23, 23.6,9),np.linspace(26.8, 27.2, 9),np.linspace(18.4, 18.9, 6),np.linspace(15.8, 16.3, 6),np.linspace(11, 11.5, 6),np.linspace(18.4, 18.7, 7)] # Octobre
# bar_tick=[np.linspace(19, 19.7, 8),np.linspace(24.5, 26.2, 11), np.linspace(26.5, 27, 6),np.linspace(19.5, 20.3, 9),np.linspace(17.1, 17.7, 7),np.linspace(10.8, 11.1, 7),np.linspace(18.7, 19, 7)] # septembre
for izone in range(len(zone)):
    print(izone, zone[izone])


    plt.rcParams["figure.figsize"] = [5, 4.50]
    plt.rcParams["figure.autolayout"] = True
    fig2, axs = plt.subplots(nrows = 1,
                                    ncols = 1, 
                                    sharex = False, 
                                    figsize=(13, 7),dpi=150) # 55,8
    df_ech = pd.read_csv(dir_lola+str(expe[0])+'_'+str(zone[izone])+'_ech_'+str(period)+'_prev-mer'+'.csv')
    df_ech.head()
    df_ech.dtypes
    #df_ech.columns.values GET name of columns
    df_ech=df_ech.loc[((df_ech['expe']==expe[0]) )][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar','temp_exp_mean','temp_merc_mean']]
    # df_ech=df_ech[5915:12239] # octobre
    df_ech=df_ech[0:5915] # septembre
    df_ech=df_ech[0:204] # septembre
    
    
    # DROP first raw
    # df_ech.drop(index=df_ech.index[0], axis=0, inplace=True)
    for niv in name:
        ech2 = np.unique(df_ech['ech'])  ## fct get only one time repeat variables
        N_level = name
        
        ## CREATE DF fill up with naans 3 cols (17 raws)
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'temp_exp_mean': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res_cmo = pd.DataFrame(data)
        
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'temp_merc_mean': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res_mer = pd.DataFrame(data)
    
        # MEAN each due date
        i = 0
        for e in ech2:
            for N in N_level:
                cmo_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['var']==N)].mean()['temp_exp_mean'])
                # print(biais_mean)
                df_res_cmo.iloc[i] = np.array([N,e,round(cmo_mean,5)])
                i = i + 1
                
                
        # MEAN each due date
        i = 0
        for e in ech2:
            for N in N_level:
                temp_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['var']==N)].mean()['temp_merc_mean'])
                # print(biais_mean)
                df_res_mer.iloc[i] = np.array([N,e,round(temp_mean,5)])
                i = i + 1


        # FIG
        df_res_mer["temp_merc_mean"]=df_res_mer["temp_merc_mean"].astype(float)
       
        N=1
        ## PLOT CONTOUR
        ax=axs
        ## DATE - XTICK
        ech=len(ech2)
        L=len(ech2)
        lst = list(np.arange(1,L+1))
        xx=lst
        list_ech_mer=[ 0,  6, 12, 18, 0,  6, 12, 18, 0,  6, 12, 18, 0,  6, 12, 18, 0]
        xind = xx[0:len(xx):N] # Pas des dates plot
        xlabels = list_ech_mer[0:len(list_ech_mer):N]

        x = lst
        y = level[0:12]
        X, Y = np.meshgrid(x, y)
        temp=[]
        temp=df_res_mer.temp_merc_mean
        temp=np.array(temp)
        
        ## Reshape bias en 2D
        reshaped_temp = temp.reshape(ech,12).T
        Z = reshaped_temp

        ## PLOT PROPRE
        ax.set_title('Zone: '+str(zone_name[izone]),fontsize=30,y=1.06)
        ax.set_xlabel('Valid time (UTC)',fontsize=30)
        ax.set_xticks(lst,
                labels=xlabels)
        axs.set_ylabel('Depth (m)',fontsize=30)
        
        # cbar_ax = fig.add_axes([0.3, -0.17, 0.5, 0.08])# [left,bottom,width,height]
        cbar_ax = fig2.add_axes([1.05, 0.095, 0.02, 0.78])
        CS1=ax.contourf(X, Y, Z,cmap='rainbow',linewidths=1)
        #CS1=ax.contourf(X, Y, Z,cmap='rainbow', levels=bar_tick[izone],linewidths=1)
        clb=fig2.colorbar(CS1, spacing='uniform',  cax=cbar_ax, orientation='vertical')
        cbar_ax.tick_params(labelsize=26)
        ax.tick_params(axis='y', labelsize=26)
        ax.tick_params(axis='x', labelsize=26)
        
        clb.ax.set_title('Temperature (°C)',fontsize=28, y=1.05)
       
        # plt.legend()
        fig2.colorbar(CS1, cax=cbar_ax, orientation='vertical')
        fig2.suptitle(str(mois)+' - averaged temperature in Mercator - '+str(expe), horizontalalignment = 'center',fontsize=28, y=1.02)
    print('min',min(temp))
    print('max',max(temp))
    print('averaged',np.mean(temp))

    figname2 = dir_fig2+str(expe[0])+'_xpmercator_column_ech_HN_'+str(zone[izone])+'_'+str(mois)+'.png'
    fig2.savefig(figname2,dpi=150, format='png',bbox_inches='tight')
    figname2 = dir_fig2+str(expe[0])+'_xpmercator_column_ech_HN_'+str(zone[izone])+'_'+str(mois)+'.pdf'
    fig2.savefig(figname2,dpi=150, format='pdf',bbox_inches='tight')

