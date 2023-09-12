#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 11:16:39 2023

@author: ormieresl
"""
from matplotlib import font_manager as fm, rcParams
import numpy.ma as ma
import os
import pandas as pd
import vortex
from datetime import datetime, timedelta
import epygram
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import usevortex
import common.util.usepygram
import olive
import common
from vortex import toolbox

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
color_GKCJ='dimgrey'
color_GO50='aqua'


## Definition color Pale
color_GL55='steelblue' ## GL55 '#648FFF'
color_GKOR='hotpink'  ## GKOR
color_GKPH= 'mediumpurple' ## GKPH '#DC267F'
color_GM6E='paleturquoise'
color_GMOT='peachpuff'
color_GN84='palegreen'#'r'
color_GO4A='rosybrown'
color_GNSR='palegreen'
color_GNYG='yellow'
color_GOJQ='orange'
color_GOJQ='peachpuff'
color_GOJQ='midnightblue'
color_GOJQ='rosybrown'
color_GO4A='steelblue'


color_GO50='darkgreen'
#import bronx
epygram.init_env()

dir_im = '/cnrm/recyf/Data/users/ormieresl/plot_param_cyclones/'
dir_df = '/home/ormieresl/Routines/Cyclones_caroles_csv/'

dir_obs= '/cnrm/recyf/Data/users/ormieresl/cyclones_nuissier/'



nom_cyclone='FIONA'
date='20220922'
colnames = ['ibbassin', 'member', 'list_ech', 'lat', 'lon', 'vmax850', 'vmax925', 'vmax10m', 'vmax50', 'vmax100m', 'vmax500m', 'vmaxgust1h', 'vmaxgust3h', 'pmin', 'zrmw34NE', 'zrmw34NW', 'zrmw34SW', 'zrmw34SE', 'zrmw64NE', 'zrmw64NW', 'zrmw64SW', 'zrmw64SE', 'sst', 'flat', 'fsen']
cyclone = 'FIONA_20220922_00.csv'

colnames_obs=['bassin','tcname','date','ventobs','pmerobs']
df_GN84 = pd.read_csv(dir_df+'GN84/glob025/'+str(cyclone), names=colnames, header=None)


df_GO4A = pd.read_csv(dir_df+'GO4A/glob025/'+str(cyclone), names=colnames, header=None)


df_GKCJ = pd.read_csv(dir_df+'GKCJ/glob025/'+str(cyclone), names=colnames, header=None)


df_GO50 = pd.read_csv(dir_df+'GO50/glob025/'+str(cyclone), names=colnames, header=None)

df_obs = pd.read_csv(dir_obs+nom_cyclone+'.csv', names=colnames_obs, header=0)

figure1, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,22],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,22],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,22],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,22],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('SST (°K)',fontsize=12)


ax2 = ax1.twinx()
ax2.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,13],linewidth=1.5,zorder=1, label = "Xp.finale",color=color_GN84)
ax2.plot(df_GO4A.iloc[:,2],df_GO4A.iloc[:,13],linewidth=1.5,zorder=1, label = "Xp.finale.buyos",color=color_GO4A)
ax2.plot(df_GO50.iloc[:,2],df_GO50.iloc[:,13],linewidth=1.5,zorder=1, label = "OML.off",color=color_GO50) 
ax2.plot(df_GKCJ.iloc[:,2],df_GKCJ.iloc[:,13],linewidth=1.5,zorder=1, label = "REF.noOML",color=color_GKCJ) 
ax2.set_ylabel('Pmin (hPa)',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
ax1.grid(alpha=0.4)
plt.title(nom_cyclone+' '+date)
plt.legend()


figure2, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,22],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,22],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,22],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,22],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('SST (°K)',fontsize=12)


ax2 = ax1.twinx()
ax2.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,7],linewidth=1.5,zorder=1, label = "Xp.finale",color=color_GN84)
ax2.plot(df_GO4A.iloc[:,2],df_GO4A.iloc[:,7],linewidth=1.5,zorder=1, label = "Xp.finale.buyos",color=color_GO4A)
ax2.plot(df_GO50.iloc[:,2],df_GO50.iloc[:,7],linewidth=1.5,zorder=1, label = "OML.off",color=color_GO50) 
ax2.plot(df_GKCJ.iloc[:,2],df_GKCJ.iloc[:,7],linewidth=1.5,zorder=1, label = "REF.noOML",color=color_GKCJ) 
ax2.set_ylabel('Vmax10',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
ax1.grid(alpha=0.4)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname2 = dir_im+nom_cyclone+'_'+date+'_SST_vmax10.png'
figure2.savefig(figname2,dpi=150, format='png',bbox_inches=None, pad_inches=0.1)

# figure2, ax1 = plt.subplots()
# ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,23],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
# ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,23],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
# ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,23],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
# ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,23],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
# ax1.set_xlabel('Forcast term',fontsize=12)
# ax1.set_ylabel('Lat',fontsize=12)

# ax2 = ax1.twinx()
# ax2.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,24],linewidth=1.5,zorder=1, label = "Xp.finale",color=color_GN84)
# ax2.plot(df_GO4A.iloc[:,2],df_GO4A.iloc[:,24],linewidth=1.5,zorder=1, label = "Xp.finale.buyos",color=color_GO4A)
# ax2.plot(df_GO50.iloc[:,2],df_GO50.iloc[:,24],linewidth=1.5,zorder=1, label = "OML.off",color=color_GO50) 
# ax2.plot(df_GKCJ.iloc[:,2],df_GKCJ.iloc[:,24],linewidth=1.5,zorder=1, label = "REF.noOML",color=color_GKCJ) 
# ax2.set_ylabel('Sens ',fontsize=12)
# xticks=np.arange(0,108,6)
# ax1.set_xticks(xticks)
# ax1.grid(alpha=0.4)
# plt.title(nom_cyclone+' '+date)
# plt.legend()


#SST
figure3, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,22],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,22],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,22],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,22],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('SST (°K)',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname3 = dir_im+nom_cyclone+'_'+date+'_SST.png'
figure3.savefig(figname3,dpi=150, format='png',bbox_inches=None, pad_inches=0.1)





#-----------------------------------------------------------------------------
# Pmin
figure4, ax2 = plt.subplots()
ax2.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,13],linewidth=1.5,zorder=1, label = "Xp.finale",color=color_GN84)
ax2.plot(df_GO4A.iloc[:,2],df_GO4A.iloc[:,13],linewidth=1.5,zorder=1, label = "Xp.finale.buyos",color=color_GO4A)
ax2.plot(df_GO50.iloc[:,2],df_GO50.iloc[:,13],linewidth=1.5,zorder=1, label = "OML.off",color=color_GO50) 
ax2.plot(df_GKCJ.iloc[:,2],df_GKCJ.iloc[:,13],linewidth=1.5,zorder=1, label = "REF.noOML",color=color_GKCJ) 
ax2.plot(df_GKCJ.iloc[:,2],df_obs.pmerobs[14:31], linewidth=1.5,zorder=1, label = "REF.noOML",color='red', linestyle='--')

ax2.set_ylabel('Pmin (hPa) ',fontsize=12)
xticks=np.arange(0,108,6)
ax2.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname4 = dir_im+nom_cyclone+'_'+date+'_Pmin_obs.png'
figure4.savefig(figname4,dpi=157, format='png',bbox_inches=None, pad_inches=0.1)





#LAT
figure5, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,23],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,23],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,23],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,23],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('Lat',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname5 = dir_im+nom_cyclone+'_'+date+'_Lat.png'
figure5.savefig(figname5,dpi=150, format='png',bbox_inches=None, pad_inches=0.1)

#SENS
figure6, ax2 = plt.subplots()
ax2.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,24],linewidth=1.5,zorder=1, label = "Xp.finale",color=color_GN84)
ax2.plot(df_GO4A.iloc[:,2],df_GO4A.iloc[:,24],linewidth=1.5,zorder=1, label = "Xp.finale.buyos",color=color_GO4A)
ax2.plot(df_GO50.iloc[:,2],df_GO50.iloc[:,24],linewidth=1.5,zorder=1, label = "OML.off",color=color_GO50) 
ax2.plot(df_GKCJ.iloc[:,2],df_GKCJ.iloc[:,24],linewidth=1.5,zorder=1, label = "REF.noOML",color=color_GKCJ) 
ax2.set_ylabel('Sens ',fontsize=12)
xticks=np.arange(0,108,6)
ax2.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname6 = dir_im+nom_cyclone+'_'+date+'_Sens.png'
figure6.savefig(figname6,dpi=150, format='png',bbox_inches=None, pad_inches=0.1)


#Vmax850
figure7, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,5],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,5],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,5],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,5],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
ax1.plot(df_GKCJ.iloc[:,2],df_obs.ventobs[14:31], linewidth=1.5,zorder=1, label = "REF.noOML",color='red', linestyle='--')

ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('Vmax850',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname7 = dir_im+nom_cyclone+'_'+date+'_Vmax850.png'
figure7.savefig(figname7,dpi=150, format='png',bbox_inches=None, pad_inches=0.1)

#Vmax10
figure8, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,7],linewidth=1.5, label = "Xp.finale",color=color_GN84)
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,7],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A)
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,7],linewidth=1.5, label = "OML.off",color=color_GO50) 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,7],linewidth=1.5, label = "REF.noOML",color=color_GKCJ) 
ax1.plot(df_GKCJ.iloc[:,2],df_obs.ventobs[14:31], linewidth=1.5,zorder=1, label = "REF.noOML",color='red', linestyle='--')

ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('Vmax10',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname8 = dir_im+nom_cyclone+'_'+date+'_Vmax10_obs.png'
figure8.savefig(figname8,dpi=170, format='png',bbox_inches=None, pad_inches=0.1)







#Vmax500
figure9, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,10],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,10],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,10],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,10],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('Vmax100',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname9 =  dir_im+nom_cyclone+'_'+date+'_Vmax100.png'
figure9.savefig(figname9,dpi=150, format='png',bbox_inches=None, pad_inches=0.1)

#Vmaxgust1h
figure10, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,11],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,11],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,11],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,11],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('Vmaxgust1h',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname10 =  dir_im+nom_cyclone+'_'+date+'_Vmaxgust1h.png'
figure10.savefig(figname10,dpi=150, format='png',bbox_inches=None, pad_inches=0.1)

#Vmaxgust3h
figure11, ax1 = plt.subplots()
ax1.plot(df_GN84.iloc[:,2],df_GN84.iloc[:,12],linewidth=1.5, label = "Xp.finale",color=color_GN84,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO4A.iloc[:,12],linewidth=1.5, label = "Xp.finale.buyos",color=color_GO4A,linestyle='--')
ax1.plot(df_GN84.iloc[:,2],df_GO50.iloc[:,12],linewidth=1.5, label = "OML.off",color=color_GO50,linestyle='--') 
ax1.plot(df_GN84.iloc[:,2],df_GKCJ.iloc[:,12],linewidth=1.5, label = "REF.noOML",color=color_GKCJ,linestyle='--') 
ax1.set_xlabel('Forcast term',fontsize=12)
ax1.set_ylabel('Vmaxgust3h',fontsize=12)
xticks=np.arange(0,108,6)
ax1.set_xticks(xticks)
plt.grid(alpha=0.6)
plt.title(nom_cyclone+' '+date)
plt.legend()

figname11 =  dir_im+nom_cyclone+'_'+date+'_Vmaxgust3h.png'
figure11.savefig(figname11,dpi=150, format='png',bbox_inches=None, pad_inches=0.1)


