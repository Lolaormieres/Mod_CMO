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
chemin = '/home/ormieresl/Scripts/Dataframe_mer/'
ech=cterm

dir_carole='/cnrm/recyf/Work/temp/labadie/pourLola/'

dir_lola='DF_temp_expe/DF_temp_expe'

dir_lola='/home/ormieresl/Routines/DF_ech_vf_mapfactor/'

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
#Only one experience - each zone
zone_carole='glob'
expe_name='Xp.final'
zone_car=['eurat','nordat','hn20','tropiques','hs20','glob','med']
zone_name=['Eurat','North-Atlantic','North Hem.','Tropics','South Hem.','Global','Med']
# xp='GN84'
expe=['GL55','GKOR','GKPH','GN84']
expe=['GKPH']

c933_mercator='GN3C'

## FIGURE
plt.rcParams["figure.figsize"] = [5, 4.50]
plt.rcParams["figure.autolayout"] = True
fig, axs = plt.subplots(nrows = 1,
                                    ncols = 7, 
                                    sharex = False, 
                                    figsize=(55, 8),dpi=200)

# LOOP XP, mean due date
for izone, zone in enumerate(zone_car):
    print(izone, zone)

    df_ech = pd.read_csv(dir_lola+str(expe[0])+'_'+str(zone)+'_ech_septembre_octobre_prev-mer'+'.csv')
    df_ech.head()
    df_ech.dtypes
    #df_ech.columns.values GET name of columns
    df_ech=df_ech.loc[((df_ech['expe']==expe[0]) )][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
    


    # DROP first raw
    df_ech.drop(index=df_ech.index[0], axis=0, inplace=True)
    for niv in name:
        # print(niv)
        ech2 = np.unique(df_ech['ech'])  ## fct get only one time repeat variables
        N_level = name
        
        ## CREATE DF fill up with naans 3 cols (17 raws)
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'eqm': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res = pd.DataFrame(data)
    
        # MEAN each due date
        i = 0
        for e in ech2:
            for N in N_level:
                eqm_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['var']==N)].mean()['eqm'])
                # print(biais_mean)
                df_res.iloc[i] = np.array([N,e,round(eqm_mean,5)])
                i = i + 1
            
        # print('min',df_res.biais.min())
        # print('max',df_res.biais.max())
        
        # FIG
        
        
        df_res["eqm"]=df_res["eqm"].astype(float)
        # print('Resultat',df_res)
       
        N=1
        ## PLOT CONTOUR
        ax=axs[izone]
        ## DATE - XTICK
        ech=len(ech2)
        ## xLabels - Dates
        L=len(ech2)
        lst = list(np.arange(1,L+1))
        xx=lst
        list_ech=np.arange(0,102,6)
        xind = xx[0:len(xx):N] # Pas des dates plot
        xlabels = list_ech[0:len(list_ech):N]
        # xind = xx[0:len(xx):N] # Pas des dates plot
        # print(xind)
        
        x = lst
        y = level[0:12]
        X, Y = np.meshgrid(x, y)
        
        eqm=[]
        eqm=df_res.eqm
        eqm=np.array(eqm)
        
        
        ## Reshape bias en 2D
        reshaped_eqm = eqm.reshape(ech,12).T
        Z = reshaped_eqm
        
        # plt.xticks(xind,labels=xlabels,rotation=45,fontsize=10)
        ax.set_title('Zone: '+str(zone_name[izone]),fontsize=36,y=1.06)
        ax.set_xlabel('Hours',fontsize=26)
        # ax.set_xticks(lst,
        #         labels=np.arange(0,102,6))
        ax.set_xticks(lst,
                labels=xlabels)
        ax.set_ylabel('Depth (m)',fontsize=26)
    
        
        ## PLOT PROPRE

        bar_tick=np.linspace(0, 0.7, 10) # GL55
        # bar_tick=np.linspace(0, 0.9, 10) # GKOR
        bar_tick=np.linspace(0, 0.7, 15)
        
        cbar_ax = fig.add_axes([0.4, -0.17, 0.2, 0.08])
        CS1=ax.contourf(X, Y, Z,cmap='magma', levels=bar_tick,linewidths=1) ## TROPIQUES
        clb=fig.colorbar(CS1, spacing='uniform', cax=cbar_ax, orientation='horizontal')
        cbar_ax.tick_params(labelsize=22)
        ax.tick_params(axis='y', labelsize=22)
        ax.tick_params(axis='x', labelsize=22)
        clb.ax.set_title('Std (°C)',fontsize=26)
        
        # [left,bottom,width,height]
        # plt.legend()
        fig.colorbar(CS1, cax=cbar_ax, orientation='horizontal')
        # fig.suptitle('Evolution du biais de température par rapport à la Reference Mercator -'+str(expe), horizontalalignment = 'center',fontsize=15)
        fig.suptitle('Evolution of the std of temperature over the 102h forecast period againt the REF. Instananeous Mercator - '+str(expe), horizontalalignment = 'center',fontsize=36, y=1.02)
    print('min',min(eqm))
    print('max',max(eqm))
    print('averaged',np.round(np.mean(eqm),4))
    

figname = '/cnrm/recyf/Data/users/ormieresl/plot_score_water_column/'+str(expe[0])+'_eqm_column_ech.png'
fig.savefig(figname,dpi=150, format='png',bbox_inches='tight')

figname = '/cnrm/recyf/Data/users/ormieresl/plot_score_water_column/'+str(expe[0])+'_eqm_column_ech.pdf'
fig.savefig(figname,dpi=150, format='pdf',bbox_inches='tight')

#%%
#Only one zone - each experience
zone_car=['eurat','nordat','hn20','tropiques','hs20','glob']
zone_name=['Eurat','North-Atlantic','North Hem.','Tropics','South Hem.','Global']
# FOR plot
zone=['eurat']
expe=['GKOR','GL55']

## FIGURE
plt.rcParams["figure.figsize"] = [5, 4.50]
plt.rcParams["figure.autolayout"] = True
fig, axs = plt.subplots(nrows = 1,
                                    ncols = 2, 
                                    sharex = False, 
                                    figsize=(19, 8),dpi=200)
# LOOP XP, mean due date
for e in range(len(expe)):
    print(e, expe[e])

    df_ech = pd.read_csv(dir_lola+str(expe[e])+'_'+str(zone[0])+'_ech_septembre_octobre_prev-mer'+'.csv')
    df_ech.head()
    df_ech.dtypes
    df_ech=df_ech.loc[((df_ech['expe']==expe[e]) )][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
    
    # DROP first raw
    df_ech.drop(index=df_ech.index[0], axis=0, inplace=True)
    for niv in name:
        # print(niv)
        ech2 = np.unique(df_ech['ech'])  ## fct get only one time repeat variables
        N_level = name
        
        ## CREATE DF fill up with naans 3 cols (17 raws)
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'eqm': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res = pd.DataFrame(data)
    
        # MEAN each due date
        i = 0
        for h in ech2:
            for N in N_level:
                eqm_mean = float(df_ech.loc[(df_ech['ech']==h) & (df_ech['var']==N)].mean()['eqm'])
                # print(eqm_mean)
                df_res.iloc[i] = np.array([N,h,round(eqm_mean,5)])
                i = i + 1

        # FIG
        df_res["eqm"]=df_res["eqm"].astype(float)

        N=1
        ## PLOT CONTOUR
        ax=axs[e]
        ## DATE - XTICK
        ech=len(ech2)
        ## xLabels - Dates
        L=len(ech2)
        lst = list(np.arange(1,L+1))
        xx=lst
        list_ech=np.arange(0,102,6)
        xind = xx[0:len(xx):N] # Pas des dates plot
        xlabels = list_ech[0:len(list_ech):N]
        # xind = xx[0:len(xx):N] # Pas des dates plot
        # print(xind)

        
        x = lst
        y = level[0:12]
        X, Y = np.meshgrid(x, y)
        
        eqm=[]
        eqm=df_res.eqm
        eqm=np.array(eqm)
        
        ## Reshape bias en 2D
        reshaped_eqm = eqm.reshape(ech,12).T
        Z = reshaped_eqm

        ax.set_title('Expe: '+str(expe[e]),fontsize=36,y=1.06)
        ax.set_xlabel('Hours',fontsize=26)
 
        ax.set_xticks(lst,
                labels=xlabels)
        ax.set_ylabel('Depth (m)',fontsize=26)

        
    ## PLOT PROPRE
        bar_tick=np.linspace(0,0.6,7)
        # bar_tick=lev # GN84
    
        cbar_ax = fig.add_axes([0.2, -0.17, 0.6, 0.08]) # [left,bottom,width,height]    
        
        CS1=ax.contourf(X, Y, Z,cmap='magma', levels=bar_tick,linewidths=1) ## TROPIQUES
        clb=fig.colorbar(CS1, spacing='uniform', cax=cbar_ax, orientation='horizontal')
        cbar_ax.tick_params(labelsize=22)
        ax.tick_params(axis='y', labelsize=22)
        ax.tick_params(axis='x', labelsize=22)
        clb.ax.set_title('Std (°C)',fontsize=26)
        
        
       
        fig.colorbar(CS1, cax=cbar_ax, orientation='horizontal')
        # fig.suptitle('Evolution du eqm de température par rapport à la Reference Mercator -'+str(expe), horizontalalignment = 'center',fontsize=15)
        fig.suptitle('Evolution of the std of temperature over the 102h \n forecast period againt the REF. Instananeous Mercator - '+str(zone[0]), horizontalalignment = 'center',fontsize=36, y=1.1)
    # print('min',min(eqm))
    # print('max',max(eqm))
    # print('averaged', np.round(np.mean(eqm),4))
    
    
figname = '/cnrm/recyf/Data/users/ormieresl/plot_score_water_column/Scores_echeance_water_columns_'+str(zone[0])+'_'+str(expe[0])+'_'+str(expe[1])+'_std.png'
fig.savefig(figname,dpi=150, format='png',bbox_inches='tight')

 
figname = '/cnrm/recyf/Data/users/ormieresl/plot_score_water_column/Scores_echeance_water_columns_'+str(zone[0])+'_'+str(expe[0])+'_'+str(expe[1])+'_std.pdf'
fig.savefig(figname,dpi=150, format='pdf',bbox_inches='tight')