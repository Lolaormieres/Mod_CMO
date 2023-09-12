#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 14:51:19 2023

@author: ormieresl
"""


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

import pandas as pd


dir_csv='/cnrm/recyf/Data/users/ormieresl/plot_score_water_column_temp_v3/'
dir_fig='/cnrm/recyf/Data/users/ormieresl/plot_score_mercator_echeances/profil_verticaux/'


## Definition color
color_GL55='#648FFF'  ## GL55
color_GKOR='#785EF0'  ## GKOR
color_GKPH='red'  ## GKPH '#DC267F'
color_GM6E='#FE6100'  ## GM6E
color_GMOT='#009E73'  ## GMOT
color_GK7C='dimgrey' 
color_GN3C='gold'
color_GN84='palegreen'
color_GNBA='midnightblue'

color_merc='gold'

zone=['eurat','med','tropiques','nordat','hn20','hs20','glob']
zone_name=['Eurat','Mediteranean','Tropics', 'North Atlantic', 'North Hemisphere', 'South Hemisphere', 'Global' ]
# zone=['eurat']
# zone_name=['Eurat']

# ## CMO
# for izone in range(len(zone)):
#     print(izone, zone[izone])
        
#     df_cmo = pd.read_csv(dir_csv+str(zone[izone])+'_CMO_september_octobre.csv', index_col=0)
#     df_cmo.head()
#     df_cmo.dtypes
    
    
#     # plt.rcParams["figure.figsize"] = [5, 4.50]
#     # plt.rcParams["figure.autolayout"] = True
#     fig, axs = plt.subplots(nrows = 1,
#                                 ncols = 1, 
#                                 sharex = False, figsize=(12,9),dpi=150)
    
    
#     df_cmo=df_cmo.T
#     x=np.arange(0,102,6)
#     heure=[0,6,12,18,0,6,12,18,0,6,12,18,0,6,12,18,0]
    
#     new_column_names = ['temp']
#     df_cmo.columns = new_column_names
    
#     df_cmo.insert(0, "ech", x, True)
#     df_cmo.insert(1, "heure", heure, True)
#     df_cmo.insert(3, "DTR", np.nan, True)
    
    
    
#     axs.plot(df_cmo.ech,df_cmo.temp,label='Xp.OML',color=color_GN84,linewidth=3)
    
#     plt.title(zone[1])
    
#     test1 = df_cmo.temp[3]-df_cmo.temp[1]
#     test2 = df_cmo.temp[7]-df_cmo.temp[5]
#     test3 = df_cmo.temp[11]-df_cmo.temp[9]
#     test4 = df_cmo.temp[15]-df_cmo.temp[13]
    
    
#     df_cmo.DTR[3]=test1
#     df_cmo.DTR[7]=test2
#     df_cmo.DTR[11]=test3
#     df_cmo.DTR[15]=test4
    
#     mean_DTR_cmo=(test1+test2+test3+test4)/4
#     print('CMOOOOOOOOOOOOOOOOOOOOOO DTR:', mean_DTR_cmo)
    
    
    
    
    
#     ## MERCATOR
#     df_merc = pd.readdf_cmo = pd.read_csv(dir_csv+str(zone[izone])+'_MERC_september_octobre.csv', index_col=0,nrows=1)
#     df_merc.head()
#     df_merc.dtypes
    
    
#     # # plt.rcParams["figure.figsize"] = [5, 4.50]
#     # # plt.rcParams["figure.autolayout"] = True
#     # fig, axs = plt.subplots(nrows = 1,
#     #                             ncols = 1, 
#     #                             sharex = False, figsize=(10,9),dpi=150)
    
    
#     df_merc=df_merc.T
#     x=np.arange(0,102,6)
#     heure=[0,6,12,18,0,6,12,18,0,6,12,18,0,6,12,18,0]
    
#     new_column_names = ['temp']
#     df_merc.columns = new_column_names
    
#     df_merc.insert(0, "ech", x, True)
#     df_merc.insert(1, "heure", heure, True)
#     df_merc.insert(3, "DTR", np.nan, True)
    
    
    
#     axs.plot(df_merc.ech,df_merc.temp,label='Mercator.instantaneous',color=color_merc,linewidth=3)
    
#     plt.title(zone[1]+'Merc')
    
#     test1 = df_merc.temp[3]-df_merc.temp[1]
#     test2 = df_merc.temp[7]-df_merc.temp[5]
#     test3 = df_merc.temp[11]-df_merc.temp[9]
#     test4 = df_merc.temp[15]-df_merc.temp[13]
    
    
#     df_merc.DTR[3]=test1
#     df_merc.DTR[7]=test2
#     df_merc.DTR[11]=test3
#     df_merc.DTR[15]=test4
    
#     mean_DTR_merc=(test1+test2+test3+test4)/4
#     print('MERCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC DTR:', mean_DTR_merc)
#     axs.set_xticks(df_cmo.ech[np.arange(17)],
#                         labels=np.arange(0,102,6))
         
    
#     axs.set_ylabel('SST (°C)',fontsize=23)
#     axs.set_xlabel('Forecast range (hours)',fontsize=24)
#     axs.set_title(str(zone_name[izone]), fontsize = 33, y= 1.01)
      
#     axs.tick_params(axis='y', labelsize=20)
#     axs.tick_params(axis='x', labelsize=20) 
#     # ax.yaxis.set_label_coords(-0.1, .5) # place of labels
#     # ax.xaxis.set_label_coords(0.5, -0.15)
#     axs.tick_params(width=2, length=4)
    
#     axs.grid(alpha=0.8)
    
#     # plt.show()
#     axs.legend(fontsize=21)
    
    
    
    # figname = dir_fig+str(zone[izone]+'_Merc_vs_CMO.png')
    # fig.savefig(figname,dpi=150, format='png',bbox_inches='tight')
    # figname = dir_fig+str(zone[izone]+'_Merc_vs_CMO.pdf')
    # fig.savefig(figname,dpi=150, format='pdf',bbox_inches='tight')
    
ech='90'


    ### PROFIL
for izone in range(len(zone)):
    print(izone, zone[izone])
    
    df_cmo = pd.read_csv(dir_csv+str(zone[izone])+'_CMO_september_octobre.csv', index_col=0)
    df_cmo.head()
    df_cmo.dtypes
    
    df_merc = pd.readdf_cmo = pd.read_csv(dir_csv+str(zone[izone])+'_MERC_september_octobre.csv', index_col=0)
    df_merc.head()
    df_merc.dtypes
    #FIG
    
    
    
    
    ech2 = [ 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96] 
    new_column_names2 = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q']
    
 
        # print()
        
    # for index, value in enumerate(df_cmo[[f'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q']]):
    #     print(f'{index}: {value}')
        
    df_cmo.columns = new_column_names2
    df_merc.columns = new_column_names2
        
    # for lettre in df_cmo[[f'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q']]:
    #     print(lettre)
    
  
    for ilettre, lettre in enumerate(new_column_names2):
        print(ilettre,lettre)
        

        # , new_column_names2[lettre])
        
        # heure= new_column_names2[lettre]
        lev=np.arange(0,12,1)
        df_cmo.insert(0, "lev", lev, True)
        df_merc.insert(0, "lev", lev, True)
        
        fig, axs = plt.subplots(nrows = 1,
                            ncols = 1, 
                            sharex = False, figsize=(7,9),dpi=150)  
        axs.plot(df_cmo[str(lettre)],df_cmo.lev,label='Xp.OML',color=color_GN84,linewidth=3)
        axs.plot(df_merc[str(lettre)],df_merc.lev,label='Mercator.instantaneous',color=color_merc,linewidth=3)
        axs.legend(fontsize=21)
        
        axs.set_ylabel('Level',fontsize=23)
        axs.set_xlabel('Temperature',fontsize=24)
        axs.set_title(str(zone_name[izone])+', P'+str(ech2[ilettre]), fontsize = 31, y= 1.1)
          
        axs.tick_params(axis='y', labelsize=20)
        axs.tick_params(axis='x', labelsize=20) 
        axs.tick_params(width=2, length=4)
        
        # axs.grid(alpha=0.8)
    
        # ax.invert_xaxis()
        axs.invert_yaxis()
        axs.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
        # axs.legend(fontsize=21)
        # plt.legend()
        # plt.show()
        
        
        
        
        figname = dir_fig+str(ech)+'_'+str(zone[izone]+'_Merc_vs_CMO_profil_verticaux.png')
        fig.savefig(figname,dpi=150, format='png',bbox_inches='tight')
        figname = dir_fig+str(ech)+'_'+str(zone[izone]+'_Merc_vs_CMO_profil_verticaux.pdf')
        fig.savefig(figname,dpi=150, format='pdf',bbox_inches='tight')
        
        # plt.cla()