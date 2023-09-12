#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 17:59:55 2023

@author: ormieresl
"""
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

# Path to stock results
pathres='/home/ormieresl/Routines/DF_bouees_vf_flag/'
dir_fig='/cnrm/recyf/Data/users/ormieresl/plot_evol_sst/'
# experience to evaluate
# expe=['GN84','GO4A']
expe=['GKCJ']


# reference
expe_mer=['GO48']
compar='prev-mer'+expe_mer[0]

# Levels numbers
nb_level=12

# Start date
date_init = datetime(2022, 9,2, 0)
deltat = timedelta(days=1)


# Days numbers
nb=30
nbDate=nb # buyos dates
period = timedelta(days=nb)

# échéances
ech_fin=102
ech_ini=0
ech_step=6

# cterm_for='96'
zone =  'nordat'
flag='flag'
# lecture des variables
name_sst=['SFX.SST']
name=['SFX.TEMP_OC'+str(i) for i in range(1,nb_level+1)]
name_sal=['SFX.SALT_OC'+str(i) for i in range(1,nb_level+1)]
name_mer=['SURF.THETAO'+str(i) for i in range(1,nb_level+1)]
name_sal_mer=['SURF.SALINO'+str(i) for i in range(1,nb_level+1)]
name_u=['SFX.UCUR_OC'+str(i) for i in range(1,nb_level+1)]
name_v=['SFX.VCUR_OC'+str(i) for i in range(1,nb_level+1)]

## Lecture des niveaux de la CMO
filename = 'level2.txt'
levelfic = np.loadtxt(filename, delimiter=',', dtype=float)
level = levelfic[0:nb_level]

## Definition color
color_GL55='#648FFF'  ## GL55
color_GKOR='#785EF0'  ## GKOR
color_GKPH='red'  ## GKPH '#DC267F'
color_GM6E='#FE6100'  ## GM6E
color_GMOT='#009E73'  ## GMOT
color_GK7C='black' 
color_GKCJ='dimgrey' 
color_GN3C='gold'
color_GN84='red'
color_GNBA='lightseagreen'
color_GO4A='orange'
color_GNSR='lightseagreen'


########################################################
#definition des champs a extraire
########################################################
## ARP - assimilation
cfields = ''
cblock = 'surfan'
ckind = 'analysis'
cut = 'assim'
cmodel='surfex'
cterm='0'
formats='fa'
fill_cmo='surf'
# cblock='forecast'
# ckind='historic'

## ARP - production - prevision
cblock_for = 'forecast'
ckind_for = 'historic'
cmodel_for = 'surfex'
# cterm_for='96'
cterm_for='0'
cut_for = 'production'
formats_for='fa'
fill_for=''

## Mercator 
cblock_oc = 'c933'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
#cfields_oc = 'sst'
cmodel_oc = 'arpege'
cut_oc = 'assim'

geometry = 'global1798'

# =============================================================================
# ## Lectur ficher fa
# =============================================================================
#LECTUR FILES
def get_file(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,
          'local':'tmpbouees_double.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpbouees_double.fic',openmode='r')
    os.remove('tmpbouees_double.fic')
    return r


#%%
# =============================================================================
# BOUEES, andepar et fgdepar
# =============================================================================
hour='0h'

fig1, ax = plt.subplots(2, 1, figsize=(15, 10),dpi=200)
# fig, ax = plt.subplots(2, 1)
fig1.tight_layout(pad=28) 

colnames=['date_valid', 'date_xp', 'expe', 'level','ech', 'obs_val', 'andepar_r', 'fgdepar_r', 'Nb_bouees_tot', 'Nb_bouees_an', 'Nb_bouees_fg','Nb_diff']
#either specify names, or inder=False
df_bouees_tot = pd.read_csv('DF_bouees_vf/'+str(zone)+'_'+str(hour)+'_andepar_fgdepar_datumflag.csv',names=colnames,header=0) #hearder=0 delete first row with col names and then get column names
# Remove column name 'A'
# df_bouees_tot.head()
# df_bouees_tot.dtypes
# df.drop(['A'], axis=1)

# df_bouees_tot.drop(['Unnamed']) # delete first row

df_GKCJ_tot=df_bouees_tot.loc[((df_bouees_tot['expe']=='GKCJ'))]

df_GN84_tot=df_bouees_tot.loc[((df_bouees_tot['expe']=='GN84'))]

df_GO4A_tot=df_bouees_tot.loc[((df_bouees_tot['expe']=='GO4A'))]


N=2
L=len(df_GN84_tot.date_valid)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84_tot.date_valid[0:len(df_GN84_tot.date_valid):N]  
## ARRANGE FIGURE
## BIAIS


ax[0].set_xticks([])
ax[0].set_ylabel('-andepar (°C)',fontsize=16)
ax[0].yaxis.set_tick_params(labelsize=16)


plt.suptitle('Buyos innovation ' +str(flag) +'\n'+str(zone),fontsize=20, y=0.95)
ax[0].plot(df_GN84_tot.date_valid,df_GN84_tot.andepar_r,label='Xp.finale',color=color_GN84,linewidth=3)
ax[0].plot(df_GO4A_tot.date_valid,df_GO4A_tot.andepar_r,label='GO4A',color=color_GO4A,linewidth=3)
ax[0].plot(df_GKCJ_tot.date_valid,df_GKCJ_tot.andepar_r,label='GKCJ',color=color_GKCJ,linewidth=3)
ax[0].axhline(y=0, color='black', linestyle='--', linewidth=2)
# ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)
ax2 = ax[0].twinx()
ax2.plot(df_GN84_tot.date_valid,df_GN84_tot.Nb_bouees_an,color='lightblue',linestyle=':',linewidth=3)
ax2.set_ylabel('Nb buyos', fontsize=16, color='lightblue')
ax2.yaxis.set_tick_params(labelsize=12)
ax2.tick_params(axis='y', colors='lightblue')


## STD
xlabels = df_GN84_tot.date_valid[0:len(df_GN84_tot.date_valid):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=16)
# ax[1].set_xlabel('Date',fontsize=16)
ax[1].set_ylabel('-fgdepar (°C)',fontsize=16)
ax[1].yaxis.set_tick_params(labelsize=16)

ax[1].plot(df_GN84_tot.date_valid,df_GN84_tot.fgdepar_r,label='Xp.finale',color=color_GN84,linewidth=3)
ax[1].plot(df_GO4A_tot.date_valid,df_GO4A_tot.fgdepar_r,label='Xp.finale',color=color_GO4A,linewidth=3)
ax[1].plot(df_GKCJ_tot.date_valid,df_GKCJ_tot.fgdepar_r,label='Xp.finale',color=color_GKCJ,linewidth=3)
ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.finale','SST.cst.mer'],fontsize=14, loc='upper right')
ax[1].legend(['Xp.final','GO4A','GKCJ.double'],fontsize=14, loc='best')

ax2 = ax[1].twinx()
ax2.plot(df_GN84_tot.date_valid,df_GN84_tot.Nb_bouees_fg,color='lightblue',linestyle=':',linewidth=3)
ax2.set_ylabel('Nb buyos', fontsize=16, color='lightblue')
ax2.yaxis.set_tick_params(labelsize=12)
ax2.tick_params(axis='y', colors='lightblue')

plt.show()

figname = dir_fig+str(zone)+'_andepar_r'+'_'+str(hour)+str(flag)+'.png'
fig1.savefig(figname,dpi=300, format='png',bbox_inches='tight')
figname = dir_fig+str(zone)+'andepar_r'+'_'+str(hour)+str(flag)+'.pdf'
fig1.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')



#%%
Dir='DF_bouees_ech_vf_new/'
# =============================================================================
# FORECAST
# =============================================================================
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
# N=2


fig, ax = plt.subplots(2, 1, figsize=(12, 9),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig.tight_layout(pad=28) 

#WITH DF ONE ECH
cterm_for='6'
df_bouees_double= pd.read_csv(pathres+str(zone)+'_P'+str(cterm_for)+'_GKCJ_bouees_'+str(flag)+'.csv')
df_bouees_double.head()
df_bouees_double.dtypes
df_GKCJ=df_bouees_double.loc[((df_bouees_double['expe']=='GKCJ'))][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees_tot','Nb_bouees_an','diff_Nb']]

# cterm_for='96'
df_bouees= pd.read_csv(pathres+str(zone)+'_P'+str(cterm_for)+'_bouees_flag.csv')
df_bouees.head()
df_bouees.dtypes
# df_GL55.drop(index=df_GL55.index[0], axis=0, inplace=True)

df_GO4A=df_bouees.loc[((df_bouees['expe']=='GO4A'))][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees_tot','Nb_bouees_an','diff_Nb']]
df_GN84=df_bouees.loc[((df_bouees['expe']=='GN84'))][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees_tot','Nb_bouees_an','diff_Nb']]
#WITH DF ALL ECH
# df_bouees_GO4A = pd.read_csv(Dir+'GO4A_'+str(zone)+'_prevision_bouees_flag_ech.csv')
# df_GO4A= df_bouees_GO4A.loc[((df_bouees_GO4A['ech']==6))]

# df_bouees_GN84 = pd.read_csv(Dir+'GN84_'+str(zone)+'_prevision_bouees_flag_ech.csv')
# df_GN84 = df_bouees_GN84.loc[((df_bouees_GN84['ech']==6))]




N=2
L=len(df_GN84.date_valid)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84.date_valid[0:len(df_GN84.date_valid):N]  
## ARRANGE FIGURE
## BIAIS


ax[0].set_xticks([])
ax[0].set_ylabel('Biais (°C)',fontsize=16)
ax[0].yaxis.set_tick_params(labelsize=16)


plt.suptitle(str(cterm_for)+'h Forecasted SST CMO 1D - SST Buyos (obs value), \n '+str(zone),fontsize=18, y=0.95)
ax[0].plot(df_GN84.date_valid,df_GN84.biais_obs,label='Xp.finale',color=color_GN84,linewidth=3)
ax[0].plot(df_GO4A.date_valid,df_GO4A.biais_obs,label='GO4A',color=color_GO4A,linewidth=3)
ax[0].plot(df_GKCJ.date_valid,df_GKCJ.biais_obs,label='GKCJ',color=color_GKCJ,linewidth=3)
ax[0].axhline(y=0, color='black', linestyle='--', linewidth=2)
# ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)
## Change line opacity, alpha=0.5
ax2 = ax[0].twinx()
# ax2.plot(df_GN84.date_valid,df_GO4A.Nb_bouees_tot,color='lightgrey',linestyle=':',linewidth=3,alpha=0.5)
ax2.plot(df_GN84.date_valid,df_GO4A.Nb_bouees_an,color=color_GO4A,linestyle=':',linewidth=3,alpha=0.6)
ax2.plot(df_GN84.date_valid,df_GN84.Nb_bouees_an,color=color_GN84,linestyle=':',linewidth=3,alpha=0.6)
ax2.plot(df_GKCJ.date_valid,df_GKCJ.Nb_bouees_an,color=color_GKCJ,linestyle=':',linewidth=3,alpha=0.6)
ax2.set_ylabel('Nb buyos', fontsize=16, color='lightgrey', alpha=0.7)
ax2.yaxis.set_tick_params(labelsize=12)
ax2.tick_params(axis='y', colors='lightgrey')


## STD
xlabels = df_GN84.date_valid[0:len(df_GN84.date_valid):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=16)
# ax[1].set_xlabel('Date',fontsize=16)
ax[1].set_ylabel('Std (°C)',fontsize=16)
ax[1].yaxis.set_tick_params(labelsize=16)

ax[1].plot(df_GN84.date_valid,df_GN84.eqm,label='Xp.finale',color=color_GN84,linewidth=3)
ax[1].plot(df_GO4A.date_valid,df_GO4A.eqm,label='Xp.finale',color=color_GO4A,linewidth=3)
ax[1].plot(df_GKCJ.date_valid,df_GKCJ.eqm,label='Xp.finale',color=color_GKCJ,linewidth=3)
ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.finale','SST.cst.mer'],fontsize=14, loc='upper right')
ax[1].legend(['Xp.final','GO4A','GKCJ.double'],fontsize=14, loc='best')
plt.show()

figname = dir_fig+'SST_forecast-obsvalue'+'_'+str(cterm_for)+'_'+str(zone)+'_'+str(flag)+'.png'
fig.savefig(figname,dpi=300, format='png',bbox_inches='tight')
figname = dir_fig+'SST_forecast-obsvalue'+'_'+str(cterm_for)+'_'+str(zone)+'_'+str(flag)+'.pdf'
fig.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')


#%% 
# =============================================================================
# ANALYSE
# =============================================================================

hour=0

## FIGURE POSTER
# %matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
# N=2


fig2, ax = plt.subplots(2, 1, figsize=(15, 10),dpi=200)
# fig, ax = plt.subplots(2, 1)
fig2.tight_layout(pad=28) 


df_bouees_an_double = pd.read_csv('DF_bouees_vf/'+str(zone)+'_'+str(cterm)+'_analyse_'+str(hour)+'h-bouees_double.csv')
df_bouees_an_double.head()
df_bouees_an_double.dtypes
df_GKCJ_an=df_bouees_an_double.loc[((df_bouees_an_double['expe']=='GKCJ'))][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees','diff_Nb']]


df_bouees_an= pd.read_csv('DF_bouees_vf/'+str(zone)+'_'+str(cterm)+'_analyse_'+str(hour)+'h-bouees.csv')
df_bouees_an.head()
df_bouees_an.dtypes
# df_GL55.drop(index=df_GL55.index[0], axis=0, inplace=True)

df_GO4A_an=df_bouees_an.loc[((df_bouees_an['expe']=='GO4A'))][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees','diff_Nb']]
df_GN84_an=df_bouees_an.loc[((df_bouees_an['expe']=='GN84'))][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees','diff_Nb']]


N=2
L=len(df_GN84_an.date_valid)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84_an.date_valid[0:len(df_GN84_an.date_valid):N]  
## ARRANGE FIGURE
## BIAIS


ax[0].set_xticks([])
ax[0].set_ylabel('Biais (°C)',fontsize=16)
ax[0].yaxis.set_tick_params(labelsize=16)

plt.suptitle(str(hour)+'h Analyse SST CMO 1D - SST Buyos (obs value), \n '+str(zone),fontsize=18, y=0.95)
ax[0].plot(df_GN84_an.date_valid,df_GN84_an.biais_obs,label='Xp.finale',color=color_GN84,linewidth=3)
ax[0].plot(df_GO4A_an.date_valid,df_GO4A_an.biais_obs,label='GO4A',color=color_GO4A,linewidth=3)
ax[0].plot(df_GKCJ_an.date_valid,df_GKCJ_an.biais_obs,label='GKCJ',color=color_GKCJ,linewidth=3)
ax[0].axhline(y=0, color='black', linestyle='--', linewidth=2)
ax2 = ax[0].twinx()
ax2.plot(df_GN84_an.date_valid,df_GO4A_an.Nb_bouees,color=color_GO4A,linestyle=':',linewidth=3)
ax2.set_ylabel('Nb buyos', fontsize=16, color='lightblue')
ax2.yaxis.set_tick_params(labelsize=12)
ax2.tick_params(axis='y', colors='lightblue')


## STD
xlabels = df_GN84_an.date_valid[0:len(df_GN84_an.date_valid):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=16)
# ax[1].set_xlabel('Date',fontsize=16)
ax[1].set_ylabel('Std (°C)',fontsize=16)
ax[1].yaxis.set_tick_params(labelsize=16)

ax[1].plot(df_GN84_an.date_valid,df_GN84_an.eqm,label='Xp.finale',color=color_GN84,linewidth=3)
ax[1].plot(df_GO4A_an.date_valid,df_GO4A_an.eqm,label='Xp.finale',color=color_GO4A,linewidth=3)
ax[1].plot(df_GKCJ_an.date_valid,df_GKCJ_an.eqm,label='Xp.finale',color=color_GKCJ,linewidth=3)
ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.finale','SST.cst.mer'],fontsize=14, loc='upper right')
ax[1].legend(['Xp.final','GO4A','GKCJ.double'],fontsize=14, loc='best')
plt.show()

figname = dir_fig+'SST_Analyse'+str(hour)+'h-obsvalue'+'_'+str(cterm_for)+'_'+str(zone)+'.png'
fig2.savefig(figname,dpi=300, format='png',bbox_inches='tight')
figname = dir_fig+'SST_Analyse'+str(hour)+'h-obsvalue'+'_'+str(cterm_for)+'_'+str(zone)+'.pdf'
fig2.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')