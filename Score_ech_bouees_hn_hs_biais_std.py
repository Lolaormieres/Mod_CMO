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
import common, olive
import common.util.usepygram
import usevortex
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import epygram
from datetime import datetime,timedelta
import vortex
import os
import pandas as pd
import numpy.ma as ma
epygram.init_env()
#Pour traiter les fichiers au format netCDF
from netCDF4 import Dataset
#Imports pour tracer des cartes
import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.mpl.ticker as cticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter




# geometry='global798c22'
geometry = 'global1798'
# chemin='/d0/images/ormieresl/'
chemin = '/home/ormieresl/Scripts/Dataframe_mer/'

Dir ='/home/ormieresl/Routines/DF_BOUEES_mask/'
# Dir_im='DF_bouees_ech_vf_new/plot_ech/'
Dir_im='/cnrm/recyf/Data/users/ormieresl/plot_ech_bouees_v3/'
# Dir_im='/cnrm/recyf/Data/users/ormieresl/plot_ech_bouees_R10x4/'


## Definition color Pale
color_GL55 = 'steelblue'  # GL55 '#648FFF'
color_GKOR = 'hotpink'  # GKOR
color_GKPH = 'mediumpurple'  # GKPH '#DC267F'
color_GM6E = 'paleturquoise'
color_GMOT = 'peachpuff'
color_GN84 = 'palegreen'  # 'r'
color_GN84 = 'limegreen'  # 'r'
color_GO4A = 'rosybrown'
color_GNSR = 'palegreen'
color_GNYG = 'yellow'
# color_GOJQ = 'orange'
# color_GOJQ = 'peachpuff'
# color_GOJQ = 'midnightblue'
color_GOJQ = 'sandybrown'
color_GOJQ = 'rosybrown'
color_GO4A='steelblue'

color_GKCJ='dimgrey'

color_expe=color_GM6E
color_mer='gold'


#%%
zones = {"nordat": [-80, 0, 20, 70], "eurat": [-35, 45, 20, 72], "tropiques": [0, 360, -
    20, 20], "hn20": [0, 360, 20, 90], "hs20": [0, 360, -90, -20],"glob": [0, 360, -90, 90]}


zone=['eurat','med','nordat', "hn20", "tropiques","hs20","glob"]
name_zone=['Eurat', 'Mediterranea','North Atlantic', 'North Hemisphere','Tropics','Southern Hemisphere','Global']
minb=[280,50,465,1200,320,380,1900]
maxb=[291,71,486,1251,341,401,1981]

expe=['GO4A','GKCJ','GOJQ']
name_exp=[ 'Xp.OML.buoys','Xp.RefnoOML','Xp.OML.buoys.R10x4']

expe=['GN84','GO4A','GKCJ']
name_exp=[ 'Xp.OML','Xp.OML.buoys','Xp.RefnoOML']


# expe=['GN84','GKCJ']
# name_exp=[ 'Xp.OML','Xp.RefnoOML']
# 
# color=[color_GO4A,color_GKCJ,color_GOJQ]
# color=[color_GN84,color_GO4A,color_GKCJ]
color=[color_GN84,color_GO4A,color_GKCJ]
# color=[color_GN84,color_GKCJ]

### Ceate Dataframe
raw_data = {'zone': [''],
                'expe': [''],
                'biais': [''],
                'eqm': [''],
                'Nb': ['']}
df2 = pd.DataFrame(raw_data, columns = ['zone','expe','biais','eqm','Nb'])

# LOOP XP, mean due date
for izone, zon in enumerate(zone):
    print(izone,zon)
    
    ## FIGURE
    # plt.rcParams["figure.figsize"] = [11, 3.50]
    plt.rcParams["figure.autolayout"] = True
    fig, axs = plt.subplots(nrows = 1,
                                    ncols = 1, 
                                    sharey = True, 
                                    figsize=(10, 5),dpi=150) #instead 25,5
    # LOOP XP, mean due date
    for isim, sim in enumerate(expe):
        print(isim,sim)
    
        # DF
        df_ech = pd.read_csv(Dir+str(sim)+'_'+str(zon)+'_prevision_bouees_flag_ech_maskcorr_septembre.csv')
        df_ech.head()
        df_ech.dtypes
    
        #df_ech.columns.values GET name of columns
        df_ech=df_ech.loc[((df_ech['expe']==sim) )]
        
        # df_ech=df_ech[254::]
    
        # DROP first raw
        df_ech.drop(index=df_ech.index[0], axis=0, inplace=True)
        ech2 = np.unique(df_ech['ech'])  ## fct get only one time repeat variables
        N_level = ['SFX.TEMP_OC1']
        
        ## CREATE DF fill up with naans 3 cols (17 raws): Biais
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'biais_obs': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res = pd.DataFrame(data)
        
        ## CREATE DF fill up with naans 3 cols (17 raws): std
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'eqm': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res_std = pd.DataFrame(data)
        
        ## CREATE DF fill for buyos numbers, up with naans 3 cols (17 raws): Nb buoys
        data2 = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
            'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
            'Nb_obs': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_buoys = pd.DataFrame(data2)
        
        
        # MEAN each due date: Biais
        i = 0
        for e in ech2:
            for N in N_level:
                biais_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['level']==-0.494025)].mean()['biais_obs'])
                # print(biais_mean)
                df_res.iloc[i] = np.array([N,e,round(biais_mean,5)])
                i = i + 1
                
        # MEAN each due date: Biais
        i = 0
        for e in ech2:
            for N in N_level:
                biais_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['level']==-0.494025)].mean()['eqm'])
                # print(biais_mean)
                df_res_std.iloc[i] = np.array([N,e,round(biais_mean,5)])
                i = i + 1
                
                
        # MEAN each due date: Buyos numbers
        i = 0
        for e in ech2:
            for N in N_level:
                Nb_bouees_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['level']==-0.494025)].mean()['Nb_bouees_an'])
                # print(biais_mean)
                df_buoys.iloc[i] = np.array([N,e,round(Nb_bouees_mean,3)])
                i = i + 1
                
        print('min of biais',df_res.biais_obs.min())
        min=float(df_res.biais_obs.min())
        print('max of biais',df_res.biais_obs.max())
        max=float(df_res.biais_obs.max())
        
        
        print('min of buoys numbers',df_buoys.Nb_obs.min())
        min=float(df_buoys.Nb_obs.min())
        print('max of buoys numbers',df_buoys.Nb_obs.max())
        max=float(df_buoys.Nb_obs.max())

        df_res["biais_obs"]=df_res["biais_obs"].astype(float)
        df_res_std["eqm"]=df_res_std["eqm"].astype(float)
        df_buoys["Nb_obs"]=df_buoys["Nb_obs"].astype(float)
     
        
        ## ax
        # ax=axs[izone]
        ax=axs # names axes for each zone
        
        ax.plot(df_res.ech,df_res.biais_obs,color=color[isim],linewidth=2.6, label=name_exp[isim])
        ax.plot(df_res_std.ech,df_res_std.eqm,color=color[isim],linewidth=2.6, linestyle='--')
        
        ax.set_xticks(df_res.ech[np.arange(17)],
                    labels=np.arange(0,102,6))
     
     
        axs.set_ylabel('°C',fontsize=18)
        ax.set_xlabel('Forecast range (hours)',fontsize=16)
        ax.set_title(str(name_zone[izone]), fontsize = 22, y= 1.05)
      
        ax.tick_params(axis='y', labelsize=16)
        ax.tick_params(axis='x', labelsize=16) 
        # ax.yaxis.set_label_coords(-0.13, .5) # place of labels
        # ax.xaxis.set_label_coords(0.5, -0.15)
        ax.tick_params(width=2, length=4)
        ax.grid(alpha=0.4)
        # ax.set_ylim([-0.2, 1.5]) #instead 1.55
    
        #df_buoys.Nb_obs.astype(int)
      
        # Right Y axis : buoys information
        
        #diff_Number of buoys
        ax2 = axs.twinx()
        ax2.plot(df_res.ech,df_buoys.Nb_obs.astype(int),color=color[isim], linestyle=':',linewidth=1.6)
        stepsize2=200
        start2=0
        end2=2200
        stepsize2=12
        start2=minb[izone]
        end2=maxb[izone]
        ax2.yaxis.set_ticks(np.arange(start2, end2, stepsize2))
        ax2.tick_params(axis='y', labelsize=15)
        ax2.yaxis.set_ticks([])
        # ax2.tick_params([])
       
        

        #set title second axis y
        ax3 = axs.twinx()
        ax3.yaxis.set_ticks([])
        ax3.set_ylabel('Buoys number',fontsize=17)
        ax3.tick_params(axis='y', labelsize=15)
        stepsize2=200
        start2=0
        end2=2200
        stepsize2=12
        start2=minb[izone]
        end2=maxb[izone]
        
        ax3.yaxis.set_ticks(np.arange(start2, end2, stepsize2))
        # ax3.yaxis.set_label_coords(1.15, 0.5)
        ax3.set_ylim([start2, end2])
    
        # ax.axhline(y=0,color='grey',linewidth=0.95)
        
        
        biais_moy=np.round(df_res.biais_obs.mean(),3)
        std_moy=np.round(df_res_std.eqm.mean(),3)
        Nb_moy=np.round(df_buoys.Nb_obs.mean())
    
        new_row = {'zone':zone[izone],'expe':expe[isim],'biais':biais_moy, 'eqm':std_moy,'Nb':Nb_moy}    
        df2.loc[len(df2)] = new_row
      
        # SUPRIMER 1er COL DF
        #df2.drop(index=df2.index[0], axis=0, inplace=True) ## DELETE FIR
        df2.to_csv(Dir_im+'forecast_sst_bouees.csv')
    
        
        # ax.legend(fontsize=14,loc='best') 
     
        # fig.savefig(Dir_im+'biais_std_bouees_allzones_mask_HN_GOJQ.png',bbox_inches="tight")
    fig.savefig(Dir_im+str(zone[izone])+'_biais_std_bouees_allzones_mask_vertpetard_2.png',bbox_inches="tight",dpi=160)
        
# #%%
# zone=["tropiques","hs20", "glob"]
# name_zone=['Tropics', 'South hemisphere', 'Global']
# # expe=['GL55','GKOR','GKPH','GN84', 'GO4A']
# # expe=['GN84','GO4A']
# # # color=[color_GL55,color_GKOR,color_GKPH,color_GN84,color_GO4A]
# # color=[color_GN84,color_GO4A]

# ## FIGURE
# fig, axs = plt.subplots(nrows = 1,
#                                     ncols = 3, 
#                                     sharey=True,
#                                     figsize=(18, 4),dpi=150)

# # LOOP XP, mean due date
# for izone, zon in enumerate(zone):
#     print(izone,zon)
#     # LOOP XP, mean due date
#     for isim, sim in enumerate(expe):
#         print(isim,sim)
    
#         # DF
#         df_ech = pd.read_csv(Dir+str(sim)+'_'+str(zon)+'_prevision_bouees_flag_ech_maskcorr_septembre.csv')
#         df_ech.head()
#         df_ech.dtypes
    
#         #df_ech.columns.values GET name of columns
#         df_ech=df_ech.loc[((df_ech['expe']==sim) )]
        
#         # df_ech=df_ech[254::]
    
#         # DROP first raw
#         df_ech.drop(index=df_ech.index[0], axis=0, inplace=True)
#         ech2 = np.unique(df_ech['ech'])  ## fct get only one time repeat variables
#         N_level = ['SFX.TEMP_OC1']
        
#         ## CREATE DF fill up with naans 3 cols (17 raws): biais
#         data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
#             'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
#             'biais_obs': np.repeat(np.nan, len(ech2)*len(N_level))}
#         df_res = pd.DataFrame(data)
        
        
               
#         ## CREATE DF fill up with naans 3 cols (17 raws): std
#         data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
#             'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
#             'eqm': np.repeat(np.nan, len(ech2)*len(N_level))}
#         df_res_std = pd.DataFrame(data)
        
#         ## CREATE DF fill for buyos numbers, up with naans 3 cols (17 raws)
#         data2 = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
#             'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
#             'Nb_obs': np.repeat(np.nan, len(ech2)*len(N_level))}
#         df_buoys = pd.DataFrame(data2)
        
        
        
        
#         # MEAN each due date: Biais
#         i = 0
#         for e in ech2:
#             for N in N_level:
#                 biais_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['level']==-0.494025)].mean()['biais_obs'])
#                 # print(biais_mean)
#                 df_res.iloc[i] = np.array([N,e,round(biais_mean,5)])
#                 i = i + 1
                
                        
#         # MEAN each due date: std
#         i = 0
#         for e in ech2:
#             for N in N_level:
#                 biais_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['level']==-0.494025)].mean()['eqm'])
#                 # print(biais_mean)
#                 df_res_std.iloc[i] = np.array([N,e,round(biais_mean,5)])
#                 i = i + 1
                
#         # MEAN each due date: Buyos numbers
#         i = 0
#         for e in ech2:
#             for N in N_level:
#                 Nb_bouees_mean = float(df_ech.loc[(df_ech['ech']==e) & (df_ech['level']==-0.494025)].mean()['Nb_bouees_an'])
#                 # print(biais_mean)
#                 df_buoys.iloc[i] = np.array([N,e,round(Nb_bouees_mean,3)])
#                 i = i + 1
                
#         print('min of biais',df_res.biais_obs.min())
#         min=float(df_res.biais_obs.min())
#         print('max of biais',df_res.biais_obs.max())
#         max=float(df_res.biais_obs.max())
        
#         print('min of buoys numbers',df_buoys.Nb_obs.min())
#         min=float(df_buoys.Nb_obs.min())
#         print('max of buoys numbers',df_buoys.Nb_obs.max())
#         max=float(df_buoys.Nb_obs.max())

#         df_res["biais_obs"]=df_res["biais_obs"].astype(float)
#         df_res_std["eqm"]=df_res_std["eqm"].astype(float)
#         df_buoys["Nb_obs"]=df_buoys["Nb_obs"].astype(float)
  
        
#         ## ax
#         ax=axs[izone] # names axes for each zone
        
#         ax.plot(df_res.ech,df_res.biais_obs,color=color[isim],linewidth=2)
#         ax.plot(df_res_std.ech,df_res_std.eqm,color=color[isim],linewidth=2, linestyle='--')
        
#         ax.set_xticks(df_res.ech[np.arange(17)],labels=np.arange(0,102,6))
#         # stepsize=0.05
#         # start=0.05
#         # end=0.40
#         # ax.yaxis.set_ticks(np.arange(start, end, stepsize))
#         axs[0].set_ylabel('°C',fontsize=12)
#         ax.set_xlabel('Forecast range (hours)',fontsize=12)
#         ax.set_title(str(name_zone[izone]), fontsize = 15, y= 1.01)

#         ax.tick_params(axis='y', labelsize=13)
#         ax.tick_params(axis='x', labelsize=12) 
#         ax.yaxis.set_label_coords(-0.14, .5) # place of labels
#         ax.xaxis.set_label_coords(0.5, -0.15)
#         ax.tick_params(width=1, length=2.5)
#         ax.grid(alpha=0.8)
        
#         ax.set_ylim([-0.2, 1.5])
 
#         #diff_Number of buoys
#         ax2 = axs[izone].twinx()
#         ax2.plot(df_res.ech,df_buoys.Nb_obs,color=color[isim], linestyle=':',linewidth=1.4)
#         stepsize2=400
#         start2=0
#         end2=2600
#         ax2.yaxis.set_ticks(np.arange(start2, end2, stepsize2))
#         ax2.tick_params(axis='y', labelsize=12)
#         ax2.yaxis.set_ticks([])
#         # ax2.tick_params([])
       
        
#         #set title second axis y
#         ax3 = axs[2].twinx()
#         ax3.yaxis.set_ticks([])
#         ax3.set_ylabel('Buoys number',fontsize=12)
#         ax3.tick_params(axis='y', labelsize=12)
#         stepsize2=400
#         start2=0
#         end2=2600
#         ax3.yaxis.set_ticks(np.arange(start2, end2, stepsize2))
#         ax3.yaxis.set_label_coords(1.2, 0.5)
        
        
          
#         # ax.axhline(y=0,color='grey',linewidth=0.85)
#         fig.tight_layout(w_pad=10)
#         # fig.suptitle('Two months mean biais OML mod. against in situ data',fontsize=40, y= 1.07)
#         # fig.legend(['Xp.finale bias','Xp.finale std','Xp.finale.buyos bias','Xp.finale.buyos std','Xp.R10x4 bias','Xp.R10x4 std', 'Ref.noOML bias', 'Ref.noOML std'],fontsize=10,bbox_to_anchor=(0.98, -0.02))
#         # fig.legend(['Xp.finale bias','Xp.finale std','Xp.finale.buyos bias','Xp.finale.buyos std','Xp.R10x4 bias','Xp.R10x4 std', 'Ref.noOML bias', 'Ref.noOML std'],fontsize=13,bbox_to_anchor=(1.14, 0.4))
#         fig.legend(['Xp.OML bias','Xp.OML std','Xp.OML.buyos bias','Xp.OML.buyos std', 'Ref.noOML bias', 'Ref.noOML std'],fontsize=13,bbox_to_anchor=(1.14, 0.4))
#         plt.legend() 
#         fig.savefig(Dir_im+'biais_std_bouees_allzones_mask_HS_glo_REF_R10x4.png',bbox_inches="tight",dpi=160)
         