#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:26:13 2023

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
import matplotlib.axes as ax
import cartopy.crs as ccrs
import cartopy
import epygram
from datetime import datetime, timedelta
import vortex
import pandas as pd
import os
import numpy.ma as ma

import pandas as pd
epygram.init_env()


## Definition color
color_GL55='royalblue' ## GL55 '#648FFF'
color_GKOR='#785EF0'  ## GKOR
color_GKPH= 'chocolate' ## GKPH '#DC267F'
color_GM6E='#FE6100'  ## GM6E
color_GMOT='limegreen'  ## GMOT '#009E73'
color_GK7C='black' 
color_GKCJ='dimgrey' 
color_OPER='r'
color_GN3C='gold'
color_GN84='red'
color_GNBA='royalblue'
color_GO4A='orange'
color_GNYG='limegreen'
color_GNSR='limegreen'


## Definition color Pale
color_GL55='steelblue' ## GL55 '#648FFF'
color_GKOR='hotpink'  ## GKOR
color_GKPH= 'mediumpurple' ## GKPH '#DC267F'
color_GM6E='paleturquoise'
color_GMOT='peachpuff'
color_GN84='r'
color_GO4A='rosybrown'
color_GNSR='palegreen'
color_GNYG='yellow'

cterm='0'
cterm_for='96'
# niv='SFX.SST'
niv='SFX.SST'


dir_df='/home/ormieresl/Routines/DF_bouees_ech_vf/'
dir_im='/cnrm/recyf/Data/users/ormieresl/plot_forecast_bouees/'

 
ech=24
zone='nordat'
#============================================================================
# ### Read DF an-mer
# =============================================================================
#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [0,360,-20,20],"hn20":[0,360,20,90],"hs20":[0,360,-90,-20],"glob":[0,360,-90,90],"med": [-3,16,30,44]}
compar='an-mer'

 
fig1, ax = plt.subplots(2, 1, figsize=(12, 9),dpi=300)
fig1.tight_layout(pad=28)  
# LOOP XP, mean due date
# for izone, zone in enumerate(zones):
#     print(izone, zone)

df_GL55 = pd.read_csv(dir_df+'GL55_'+str(zone)+'_prevision_bouees_flag_ech.csv')
df_GL55.head()
df_GL55.dtypes
df_GL55.drop(index=df_GL55.index[0], axis=0, inplace=True)
df_GL55_ech=df_GL55.loc[df_GL55.ech==ech]
    
df_GKOR = pd.read_csv(dir_df+'GKOR_'+str(zone)+'_prevision_bouees_flag_ech.csv')
df_GKOR.head()
df_GKOR.dtypes
df_GKOR.drop(index=df_GKOR.index[0], axis=0, inplace=True)
df_GKOR_ech=df_GKOR.loc[df_GKOR.ech==ech]


df_GKPH = pd.read_csv(dir_df+'GKPH_'+str(zone)+'_prevision_bouees_flag_ech.csv')
df_GKPH.head()
df_GKPH.dtypes
df_GKPH.drop(index=df_GKPH.index[0], axis=0, inplace=True)
df_GKPH_ech=df_GKPH.loc[df_GKPH.ech==ech]


df_GN84 = pd.read_csv(dir_df+'GN84_'+str(zone)+'_prevision_bouees_flag_ech.csv')
df_GN84.head()
df_GN84.dtypes
df_GN84.drop(index=df_GN84.index[0], axis=0, inplace=True)
df_GN84_ech=df_GN84.loc[df_GN84.ech==ech]

   
df_GO4A = pd.read_csv(dir_df+'GO4A_'+str(zone)+'_prevision_bouees_flag_ech.csv')
df_GO4A.head()
df_GO4A.dtypes
df_GO4A.drop(index=df_GO4A.index[0], axis=0, inplace=True)
df_GO4A_ech=df_GO4A.loc[df_GO4A.ech==ech]


# Plot dates axe x
N=2

L=len(df_GN84_ech.date_valid)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84_ech.date_valid[0:len(df_GN84_ech.date_valid):N]

## ARRANGE FIGURE
## BIAIS
plt.suptitle('Forecast SST (ech:'+str(ech)+')- Buyos flag \n'+str(zone),fontsize=26, y=0.97)
ax[0].set_xticks([])
ax[0].set_ylabel('Bias (°C)',fontsize=22)
ax[0].yaxis.set_tick_params(labelsize=22)

   

ax[0].plot(df_GL55_ech.date_valid,df_GL55_ech.biais_obs,label='L.26',color=color_GL55,linewidth=3)
ax[0].plot(df_GKOR_ech.date_valid,df_GKOR_ech.biais_obs,label='L.12',color=color_GKOR,linewidth=3)
ax[0].plot(df_GKPH_ech.date_valid,df_GKPH_ech.biais_obs,label='R.50',color=color_GKPH,linewidth=3)
ax[0].plot(df_GN84_ech.date_valid,df_GN84_ech.biais_obs,label='Xp.finale',color=color_GN84,linewidth=3)
ax[0].plot(df_GO4A_ech.date_valid,df_GO4A_ech.biais_obs,label='Xp.buoys',color=color_GO4A,linewidth=3)
# ax[0].plot(df_GKCJ.date,df_GKCJ.biais,label='Double GKCJ',color=color_GKCJ,linewidth=3)

ax[0].axhline(y=0, color='black', linestyle='--', linewidth=2)
# ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)
    
#Nb points - stretch grid
ax2 = ax[0].twinx()

# Nb=df_GN84_ech.Nb_bouees_an.asfloat()
ax2.plot(df_GN84_ech.date_valid,df_GN84_ech.Nb_bouees_an,color=color_GN84,linestyle=':',linewidth=3,alpha=0.55)
ax2.set_ylabel('Nb_points',color='steelblue', fontsize=22)
ax2.tick_params(axis="y", labelcolor='steelblue',labelsize=19)

## STD
# xlabels = df_GNSR.date[0:len(df_GNSR.date):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=18)
# ax[1].set_xlabel('Date',fontsize=22)
ax[1].set_ylabel('Standard deviation (°C)',fontsize=22)
ax[1].yaxis.set_tick_params(labelsize=22)

ax[1].plot(df_GL55_ech.date_valid,df_GL55_ech.eqm,label='GL55',color=color_GL55,linewidth=3)
ax[1].plot(df_GKOR_ech.date_valid,df_GKOR_ech.eqm,label='GKOR',color=color_GKOR,linewidth=3)
ax[1].plot(df_GKPH_ech.date_valid,df_GKPH_ech.eqm,label='GKPH',color=color_GKPH,linewidth=3)
ax[1].plot(df_GN84_ech.date_valid,df_GN84_ech.eqm,label='Xp.finale',color=color_GN84,linewidth=3)
ax[1].plot(df_GO4A_ech.date_valid,df_GO4A_ech.eqm,label='GO4A.buyos',color=color_GO4A,linewidth=3)
# ax[1].plot(df_GKCJ.date,df_GKCJ.eqm,label='Double GKCJ',color=color_GKCJ,linewidth=3)

ax[1].legend(['L.26','L.12','R.50','Xp.finale','Xp.buyos'],fontsize=19,bbox_to_anchor=(1.1, 1.05))
# ax2.set_ylim(100000,2300000)    
# ax[0].set_ylim(-0.2, 0.6)
# ax[1].set_ylim(0, 0.7)
   
plt.show()
figname = dir_im+'SST'+'_'+str(zone)+'_ech_'+str(ech)+'.png'
fig1.savefig(figname,dpi=300, format='png',bbox_inches='tight')
# figname = dir_fig+'SST_biais'+'_'+str(zone)+'_ech'+str(ech)+'.pdf'
# fig1.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')

   