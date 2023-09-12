#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 14:34:45 2023

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


# ## Definition color
# color_GL55='#648FFF'  ## GL55
# color_GKOR='#785EF0'  ## GKOR
# color_GKPH='crimson'  ## GKPH
# color_GM6E='#FE6100'  ## GM6E
# color_GMOT='limegreen'  ## GMOT
# color_OPER='r'
# color_GK7C='black' 
# color_GN3C='gold'


## Definition color
color_GL55='royalblue' ## GL55 '#648FFF'
color_GKOR='#785EF0'  ## GKOR
color_GKPH= '#FF0000' ## GKPH '#DC267F'
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
cterm='0'
cterm_for='96'
# niv='SFX.SST'
niv='SFX.SST'
niv_GK7C='SFX.SST'
niv_GKCJ='SFX.SST'
niv_GN84='SFX.TEMP_OC1'
# zone_GK7C = 'nord_atlantique'
# zone='nord-atlantique'
# zone_GN84 = 'nord-atlantique'
# nom_zone = 'Nord-Atlantique'
zone_GK7C = 'eurat'
zone='eurat'
zone_GN84 = 'eurat'
nom_zone = 'eurat'
# zone_GK7C = 'tropiques'
# zone='tropiques'
# zone_GN84 = 'tropiques'
# nom_zone = 'Tropiques'


#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [0,360,-20,20],"hn20":[0,360,20,90],"hs20":[0,360,-90,-20],"glob":[0,360,-90,90]}

zone=[]
for key in zones:
  
    print(key)
    zone.append(key)
print(zone)

#%%
# =============================================================================
# ### Read DF an-mer
# =============================================================================
# septembre_octobre_an-mer_sst_tropiquessst.csv

compar='an-mer'
df_GL55 = pd.read_csv('DF_SST_vf/DF_SST_vf/GL55_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GL55.head()
df_GL55.dtypes
df_GL55.drop(index=df_GL55.index[0], axis=0, inplace=True)

df_GN84 = pd.read_csv('DF_SST_vf/DF_SST_vf/GN84_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GN84.head()
df_GN84.dtypes
df_GN84.drop(index=df_GN84.index[0], axis=0, inplace=True)

df_GO4A = pd.read_csv('DF_SST_vf/DF_SST_vf/GO4A_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GO4A.head()
df_GO4A.dtypes
df_GO4A.drop(index=df_GO4A.index[0], axis=0, inplace=True)

# df_GNYG = pd.read_csv('DF_SST_vf/DF_SST_vf/GNYG_'+str(zone)+'_'+str(compar)+'_sst.csv')
# df_GNYG.head()
# df_GNYG.dtypes

df_GKCJ = pd.read_csv('DF_SST_vf/DF_SST_vf/GKCJ_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GKCJ.head()
df_GKCJ.dtypes
df_GKCJ.drop(index=df_GKCJ.index[0], axis=0, inplace=True)



#%%
# FIGURE POSTER
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
# N=2

fig2, ax = plt.subplots(2, 1, figsize=(12, 9),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig2.tight_layout(pad=28)  

# Plot dates axe x
N=8

L=len(df_GN84.date)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84.date[0:len(df_GN84.date):N]

## ARRANGE FIGURE
## BIAIS
ax[0].set_xticks([])
ax[0].set_ylabel('Bias (°C)',fontsize=22)
ax[0].yaxis.set_tick_params(labelsize=22)

# plt.suptitle('SST ARPEGE analysée - SST Mercator (REF) \n\n',fontsize=18)
# plt.suptitle('Analysed SST - Instantaneous Mercator SST (REF) \n zone: '+str(nom_zone),fontsize=26, y=0.96)
plt.suptitle('Analysed SST - Instantaneous Mercator SST (REF) \n'+str(zone),fontsize=26, y=0.97)

ax[0].plot(df_GL55.date,df_GL55.biais,label='L.26',color=color_GL55,linewidth=3)
# ax[0].plot(df_GKOR.date,df_GKOR.biais,label='L.12',color=color_GKOR,linewidth=3)
# ax[0].plot(df_GKPH.date,df_GKPH.biais,label='R.50',color=color_GKPH,linewidth=3)
# ax[0].plot(df_GM6E.date,df_GM6E.biais,label='no.curr',color=color_GM6E,linewidth=3)
# ax[0].plot(df_GMOT.date,df_GMOT.biais,label='no.curr.bathy',color=color_GMOT,linewidth=3)
ax[0].plot(df_GN84.date,df_GN84.biais,label='Xp.finale',color=color_GN84,linewidth=3)
ax[0].plot(df_GO4A.date,df_GO4A.biais,label='GO4A',color=color_GO4A,linewidth=3)
# ax[0].plot(df_OPER.date,df_OPER.biais,label='OPER.OSTIA',color=color_OPER)
ax[0].plot(df_GKCJ.date,df_GKCJ.biais,label='Double GKCJ',color=color_GKCJ,linewidth=3)
ax[0].axhline(y=0, color='black', linestyle='--', linewidth=2)
# ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)

## STD
xlabels = df_GN84.date[0:len(df_GN84.date):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=22)
# ax[1].set_xlabel('Date',fontsize=22)
ax[1].set_ylabel('Standard deviation (°C)',fontsize=22)
ax[1].yaxis.set_tick_params(labelsize=22)
# ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR,linewidth=3)
ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55,linewidth=3)
# ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH,linewidth=3)
# ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E,linewidth=3)
# ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT,linewidth=3)
ax[1].plot(df_GN84.date,df_GN84.eqm,label='Xp.finale',color=color_GN84,linewidth=3)
ax[1].plot(df_GO4A.date,df_GO4A.eqm,label='GO4A',color=color_GO4A,linewidth=3)
# ax[1].plot(df_OPER.date,df_OPER.eqm,label='OPER.OSTIA',color=color_OPER)
ax[1].plot(df_GKCJ.date,df_GKCJ.eqm,label='Double GKCJ',color=color_GKCJ,linewidth=3)
# ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.finale','SST.cst.mer'],fontsize=14, loc='upper right')
ax[1].legend(['L.26','Xp.finale','GO4A','oper.2024'],fontsize=19, loc='best')
plt.show()


figname = '/cnrm/recyf/Data/users/ormieresl/SST_evol_temporelle_analyse-mercator_'+'_'+str(zone)+'.png'
fig2.savefig(figname,dpi=300, format='png',bbox_inches='tight')
figname = '/cnrm/recyf/Data/users/ormieresl/SST_evol_temporelle_analyse-mercator_'+'_'+str(zone)+'.pdf'
fig2.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')



#%%
# =============================================================================
# ### Read DF prev-mer
# =============================================================================
compar='prev-mer'
cterm='0'

df_GL55 = pd.read_csv('DF_SST_vf/DF_SST_vf/GL55_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GL55.head()
df_GL55.dtypes
df_GL55.drop(index=df_GL55.index[0], axis=0, inplace=True)

df_GN84 = pd.read_csv('DF_SST_vf/DF_SST_vf/GN84_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GN84.head()
df_GN84.dtypes
df_GN84.drop(index=df_GN84.index[0], axis=0, inplace=True)

df_GO4A = pd.read_csv('DF_SST_vf/DF_SST_vf/GO4A_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GO4A.head()
df_GO4A.dtypes
df_GO4A.drop(index=df_GO4A.index[0], axis=0, inplace=True)

df_GNYG = pd.read_csv('DF_SST_vf/DF_SST_vf/GNYG_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GNYG.head()
df_GNYG.dtypes
df_GNYG.drop(index=df_GNYG.index[0], axis=0, inplace=True)

df_GKCJ = pd.read_csv('DF_SST_vf/DF_SST_vf/GKCJ_'+str(zone)+'_'+str(compar)+'_sst.csv')
df_GKCJ.head()
df_GKCJ.dtypes
df_GKCJ.drop(index=df_GKCJ.index[0], axis=0, inplace=True)

#%%
## FIGURE POSTER
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
# N=2


fig2, ax = plt.subplots(2, 1, figsize=(12, 9),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig2.tight_layout(pad=28) 

# Plot dates axe x
N=8

L=len(df_GN84.date)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84.date[0:len(df_GN84.date):N]


# ## BIAIS

# ax[0].set_xticks([])
# ax[0].set_ylabel('Biais (°C)',fontsize=16)
# ax[0].yaxis.set_tick_params(labelsize=16)
# plt.suptitle('SST prévue par le modèle CMO 1D (échéance : 96h) - SST Mercator (REF), \n '+str(nom_zone),fontsize=18, y=0.95)
# ax[0].plot(df_GL55.date,df_GL55.biais,label='L.26',color=color_GL55,linewidth=3)
# ax[0].plot(df_GKPH.date,df_GKPH.biais,label='R.50',color=color_GKPH,linewidth=3)
# ax[0].plot(df_GM6E.date,df_GM6E.biais,label='no.curr',color=color_GM6E,linewidth=3)
# ax[0].plot(df_GMOT.date,df_GMOT.biais,label='no.curr.bathy',color=color_GMOT,linewidth=3)
# ax[0].plot(df_GKOR.date,df_GKOR.biais,label='L.12',color=color_GKOR,linewidth=3)
# ax[0].plot(df_GN84.date,df_GN84.biais,label='Xp.finale',color=color_GN84,linewidth=3)
# # ax[0].plot(df_OPER.date,df_OPER.biais,label='OPER.OSTIA',color=color_OPER)
# ax[0].plot(df_GK7C.date,df_GK7C.biais,label='SST fixe',color=color_GK7C,linewidth=3)
# ax[0].axhline(y=0, color='black', linestyle='--', linewidth=2)
# # ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)

# ## STD
# xlabels = df_GL55.date[0:len(df_GL55.date):N]
# ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=16)
# ax[1].set_xlabel('Date',fontsize=16)
# ax[1].set_ylabel('Std (°C)',fontsize=16)
# ax[1].yaxis.set_tick_params(labelsize=16)
# ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR,linewidth=3)
# ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55,linewidth=3)
# ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH,linewidth=3)
# ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E,linewidth=3)
# ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT,linewidth=3)
# ax[1].plot(df_GN84.date,df_GN84.eqm,label='GN84',color=color_GN84,linewidth=3)
# # ax[1].plot(df_OPER.date,df_OPER.eqm,label='OPER.OSTIA',color=color_OPER)
# ax[1].plot(df_GK7C.date,df_GK7C.eqm,label='MERCATOR moyen',color=color_GK7C,linewidth=3)
# ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.finale','SST.cst. sans CMO'],fontsize=14, loc='upper right')
# plt.show()
# # plt.legend(loc=(1,0),fontsize=14)
# ##GK7C mercator moyen



# BIAIS
ax[0].set_xticks([])
ax[0].set_ylabel('Bias (°C)',fontsize=22)
ax[0].yaxis.set_tick_params(labelsize=22)

# plt.suptitle('SST ARPEGE analysée - SST Mercator (REF) \n\n',fontsize=18)
# plt.suptitle(' SST prévue par le modèle CMO 1D (échéance : 96h) - SST Mercator (REF), \n '+str(nom_zone),fontsize=18, y=0.95)
# plt.suptitle('96-h forecasted SST - Instantaneous Mercator SST (REF) \n zone: '+str(nom_zone),fontsize=26, y=0.96)
plt.suptitle('96-h forecasted SST - Instantaneous Mercator SST (REF) \n'+str(zone),fontsize=22, y=0.95)

            
ax[0].plot(df_GL55.date,df_GL55.biais,label='L.26',color=color_GL55,linewidth=3)
# ax[0].plot(df_GKOR.date,df_GKOR.biais,label='L.12',color=color_GKOR,linewidth=3)
# ax[0].plot(df_GKPH.date,df_GKPH.biais,label='R.50',color=color_GKPH,linewidth=3)
# ax[0].plot(df_GM6E.date,df_GM6E.biais,label='no.curr',color=color_GM6E,linewidth=3)
# ax[0].plot(df_GMOT.date,df_GMOT.biais,label='no.curr.bathy',color=color_GMOT,linewidth=3)
ax[0].plot(df_GN84.date,df_GN84.biais,label='Xp.finale',color=color_GN84,linewidth=3)
ax[0].plot(df_GO4A.date,df_GO4A.biais,label='GO4A',color=color_GO4A,linewidth=3)
# ax[0].plot(df_GNYG.date,df_GNYG.biais,label='GNYG_wind',color=color_GNYG,linewidth=3)
# ax[0].plot(df_OPER.date,df_OPER.biais,label='OPER.OSTIA',color=color_OPER)
ax[0].plot(df_GKCJ.date,df_GKCJ.biais,label='SST fixe',color=color_GKCJ,linewidth=3)
ax[0].axhline(y=0, color='black', linestyle='--', linewidth=2)
# ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)

## STD
xlabels = df_GN84.date[0:len(df_GN84.date):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=22)
ax[1].set_xlabel('Date',fontsize=22)
ax[1].set_ylabel('Standard deviation (°C)',fontsize=22)
ax[1].yaxis.set_tick_params(labelsize=22)
# ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR,linewidth=3)
ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55,linewidth=3)
# ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH,linewidth=3)
# ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E,linewidth=3)
# ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT,linewidth=3)
ax[1].plot(df_GN84.date,df_GN84.eqm,label='Xp.finale',color=color_GN84,linewidth=3)
ax[1].plot(df_GO4A.date,df_GO4A.eqm,label='GO4A',color=color_GO4A,linewidth=3)
# ax[1].plot(df_GNYG.date,df_GNYG.eqm,label='GNYG_wind',color=color_GNYG,linewidth=3)
# ax[1].plot(df_OPER.date,df_OPER.eqm,label='OPER.OSTIA',color=color_OPER)
ax[1].plot(df_GKCJ.date,df_GKCJ.eqm,label='MERCATOR moyen',color=color_GKCJ,linewidth=3)
# ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.finale','SST.cst.mer'],fontsize=14, loc='upper right')
ax[1].legend(['L.26','Xp.finale','GO4A','GKCJ'],fontsize=19, loc='upper right')
plt.show()

figname = '/cnrm/recyf/Data/users/ormieresl/SST_evol_temporelle_prevision-mercator_'+'_'+str(zone)+'.png'
fig2.savefig(figname,dpi=300, format='png',bbox_inches='tight')
figname = '/cnrm/recyf/Data/users/ormieresl/SST_evol_temporelle_prevision-mercator_'+'_'+str(zone)+'.pdf'
fig2.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')
