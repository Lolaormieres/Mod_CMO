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
# import seaborn as sns

dir_csv='/cnrm/recyf/Data/users/ormieresl/plot_score_water_column_temp_v3/'
dir_fig='/cnrm/recyf/Data/users/ormieresl/plot_score_mercator_echeances/sst/'


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

## CMO
for izone in range(len(zone)):
    print(izone, zone[izone])
        
    df_cmo_1 = pd.read_csv(dir_csv+str(zone[izone])+'_CMO_september_octobre.csv', index_col=0,nrows=1)
    df_cmo_1.head()
    df_cmo_1.dtypes
    
    
    # plt.rcParams["figure.figsize"] = [5, 4.50]
    # plt.rcParams["figure.autolayout"] = True
    fig, axs = plt.subplots(nrows = 1,
                                ncols = 1, 
                                sharex = False, figsize=(12,9),dpi=150)
    
    
    df_cmo=df_cmo_1.T
    x=np.arange(0,102,6)
    heure=[0,6,12,18,0,6,12,18,0,6,12,18,0,6,12,18,0]
    
    new_column_names = ['temp']
    df_cmo.columns = new_column_names
    
    df_cmo.insert(0, "ech", x, True)
    df_cmo.insert(1, "heure", heure, True)
    df_cmo.insert(3, "DTR", np.nan, True)
    
    
    
    axs.plot(df_cmo.ech,df_cmo.temp,label='Xp.OML',color=color_GN84,linewidth=3)
    
    plt.title(zone[1])
    
    test1 = df_cmo.temp[3]-df_cmo.temp[1]
    test2 = df_cmo.temp[7]-df_cmo.temp[5]
    test3 = df_cmo.temp[11]-df_cmo.temp[9]
    test4 = df_cmo.temp[15]-df_cmo.temp[13]
    
    
    df_cmo.DTR[3]=test1
    df_cmo.DTR[7]=test2
    df_cmo.DTR[11]=test3
    df_cmo.DTR[15]=test4
    
    mean_DTR_cmo=(test1+test2+test3+test4)/4
    print('CMOOOOOOOOOOOOOOOOOOOOOO DTR:', mean_DTR_cmo)
    
    ## MERCATOR
    df_merc = pd.readdf_cmo = pd.read_csv(dir_csv+str(zone[izone])+'_MERC_september_octobre.csv', index_col=0,nrows=1)
    df_merc.head()
    df_merc.dtypes
    
    
    # # plt.rcParams["figure.figsize"] = [5, 4.50]
    # # plt.rcParams["figure.autolayout"] = True
    # fig, axs = plt.subplots(nrows = 1,
    #                             ncols = 1, 
    #                             sharex = False, figsize=(10,9),dpi=150)
    
    
    df_merc=df_merc.T
    x=np.arange(0,102,6)
    heure=[0,6,12,18,0,6,12,18,0,6,12,18,0,6,12,18,0]
    
    new_column_names = ['temp']
    df_merc.columns = new_column_names
    
    df_merc.insert(0, "ech", x, True)
    df_merc.insert(1, "heure", heure, True)
    df_merc.insert(3, "DTR", np.nan, True)
    
    
    
    axs.plot(df_merc.ech,df_merc.temp,label='Mercator.instantaneous',color=color_merc,linewidth=3)
    
    plt.title(zone[1]+'Merc')
    
    test1 = df_merc.temp[3]-df_merc.temp[1]
    test2 = df_merc.temp[7]-df_merc.temp[5]
    test3 = df_merc.temp[11]-df_merc.temp[9]
    test4 = df_merc.temp[15]-df_merc.temp[13]
    
    
    df_merc.DTR[3]=test1
    df_merc.DTR[7]=test2
    df_merc.DTR[11]=test3
    df_merc.DTR[15]=test4
    
    mean_DTR_merc=(test1+test2+test3+test4)/4
    print('MERCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC DTR:', mean_DTR_merc)
    axs.set_xticks(df_cmo.ech[np.arange(17)],
                        labels=np.arange(0,102,6))
         
    
    axs.set_ylabel('SST (°C)',fontsize=23)
    axs.set_xlabel('Forecast range (hours)',fontsize=24)
    axs.set_title(str(zone_name[izone]), fontsize = 33, y= 1.01)
      
    axs.tick_params(axis='y', labelsize=20)
    axs.tick_params(axis='x', labelsize=20) 
    # ax.yaxis.set_label_coords(-0.1, .5) # place of labels
    # ax.xaxis.set_label_coords(0.5, -0.15)
    axs.tick_params(width=2, length=4)
    
    axs.grid(alpha=0.8)
    
    # plt.show()
    axs.legend(fontsize=21)
    
  
    
    figname = dir_fig+str(zone[izone]+'_Merc_vs_CMO.png')
    fig.savefig(figname,dpi=150, format='png',bbox_inches='tight')
    figname = dir_fig+str(zone[izone]+'_Merc_vs_CMO.pdf')
    fig.savefig(figname,dpi=150, format='pdf',bbox_inches='tight')
    
    
       #Import seaborn
  
#     df_ech = pd.read_csv(dir_carole+str(expe[0])+'_'+str(zone[izone])+'_ech_septembre_octobre_prev-mer'+'.csv')
#     df_ech.head()
#     df_ech.dtypes
#     #df_ech.columns.values GET name of columns
#     df_ech=df_ech.loc[((df_ech['expe']==expe[0]) )][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar','temp_exp_mean','temp_merc_mean']]
    
# sns.lineplot(data=df_cmo, x="ech", y="temp")
    
    #%%
    
dir_carole='/home/ormieresl/Routines/gmapfactor/'

expe=['GN84'] 
df_tot = pd.read_csv(dir_carole+str(expe[0])+'_'+str(zone[izone])+'_ech_septembre_octobre_prev-mer'+'.csv')
df_tot.head()
df_tot.dtypes
#df_ech.columns.values GET name of columns
df_tot=df_tot.loc[((df_tot['expe']==expe[0]) )][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar','temp_exp_mean','temp_merc_mean']]

# df_tot_cmo = pd.read_csv(dir_carole+str(expe[0])+'_'+str(zone[izone])+'_ech_septembre_octobre_prev-mer'+'.csv')
# df_tot_cmo.head()
# df_tot_cmo.dtypes
# #df_ech.columns.values GET name of columns
# df_tot_cmo=df_tot_cmo.loc[((df_tot_cmo['expe']==expe[0]) )][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar','temp_exp_mean','temp_merc_mean']]




# sns.relplot(data=tips, x="total_bill", y="tip", hue="day", col="time")

fig2, axs = plt.subplots(nrows = 1,
                                ncols = 1, 
                                sharex = False, figsize=(15,9),dpi=150)
    
ech = [ 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96] 
l_merc=[]
for a in ech:
    y_error_merc=np.std(df_tot.temp_merc_mean.iloc[np.where(df_tot.ech==a)])
    print(y_error_merc)
    l_merc.append(y_error_merc)

l_cmo=[]
for a in ech:
    y_error_cmo=np.std(df_tot.temp_exp_mean.iloc[np.where(df_tot.ech==a)])
    print(y_error_cmo)
    l_cmo.append(y_error_cmo)

axs.plot(df_merc.ech,df_merc.temp,label='Mercator.instantaneous',color=color_merc,linewidth=3)
axs.plot(df_cmo.ech,df_cmo.temp,label='Xp.OML',color=color_GN84,linewidth=3)


plt.errorbar(df_merc.ech, df_merc.temp, yerr=l_merc, color=color_merc)
plt.errorbar(df_cmo.ech, df_cmo.temp, yerr=l_cmo, color=color_GN84)
## ytest


import numpy as np; np.random.seed(22)
import seaborn as sns; sns.set(color_codes=True)
x = np.linspace(0, 15, 31)
data = np.sin(x) + np.random.rand(10, 31) + np.random.randn(10, 1)
ax = sns.tsplot(data=data, err_style="ci_bars")
plt.show()