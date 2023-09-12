#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:21:35 2023

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


## Lecture fichier Netcdf
exp='GN84'

## Date - exp
nb=1
date_init = datetime(2022, 8,2, 0)
date = date_init
dt = timedelta(days=1)
date_valid = date+dt
period = timedelta(days=nb)
date_end = date+period



## Date - Bouees - stage M1
# import datetime as dt

# year=2021
# month=4 # 4 pour Avril (et non 5 Mai)
# day=1
# hour=0
# Delta=24
# nbDate=120

# Date_debut = dt.datetime(year, month, day, hour)

# timedelta = dt.timedelta(0,0,0,0,0,Delta)

# toto=[Date_debut+i*timedelta for i in range(nbDate)]

# reseaux=""
# P=5

# # BOUCLE lecture des fichiers
# for date_dt in toto:
#     date = '%4.4i%2.2i%2.2i%2.2i00'%(date_dt.year,\
#     date_dt.month,\
#     date_dt.day,\
#     date_dt.hour)
    
#     date_2 = '%2.2i/%2.2i/%2.2i %2.2iUTC'%(date_dt.day,\
#     date_dt.month,\
#     date_dt.year,\
#     date_dt.hour)

dir_nc='/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(exp)+'/20220928T0000A/surfan/'
fic_nc = dir_nc +'odb-ecma.canari.surf.nc'

dir_fig='bouees/'

myfile = Dataset(fic_nc,'r')

#Pour donner des noms explicites aux différentes variables :
odb_key_val_dic = {}
for variables in myfile.variables.keys():
#instrument_int = myfile[variables].odb_name[:myfile[variables].odb_name.index('@')]
    odb_key_val_dic[variables] = myfile[variables][:]
            #print(variables,'->',instrument_int)
            
            
ind_flag = np.where((odb_key_val_dic['col_2']>20)&
(odb_key_val_dic['col_2']<70)&
(odb_key_val_dic['col_1']>-80)&
(odb_key_val_dic['col_1']<45))
        

#data2plota = odb_key_val_dic['obsvalue'][ind_flag]
data2plota = odb_key_val_dic['col_7'][ind_flag]
lon       = odb_key_val_dic['col_1'][ind_flag]
lat       = odb_key_val_dic['col_2'][ind_flag]
datamin   = 275
datamax   = 300


val2plot='sst_obs'
instrument='bouees'

# Test colobar discontinu
# Colormap 
cmap_type = 'plasma' # le type de colormap
nb_of_bins = 10 # le nbre de couleurs dans la colormap
cmap = plt.get_cmap(cmap_type, nb_of_bins)
vmin, vmax = 275, 300
label_data="temp"

# Colorbar continu
#         cmap = matplotlib.cm.plasma
projection=ccrs.PlateCarree(central_longitude=-45)
figsize=(20,10)
fig,ax = plt.subplots(nrows=1,ncols=1,figsize = figsize,
                                        subplot_kw={'projection': projection},sharex=True,sharey=True)

im = ax.scatter(lon,lat,c=data2plota, transform=ccrs.PlateCarree(),vmin = datamin,vmax=datamax,cmap=cmap,s=50.0,marker="o")
ax.coastlines()
divider = make_axes_locatable(ax)
cax     = divider.append_axes("right", size="5%", pad=0.35, axes_class=plt.Axes)
## Colorbar continu
cbar     = plt.colorbar(im,cax=cax)
plt.ylabel('SST (K)', fontsize=20)
ticklabs = cbar.ax.get_yticklabels()
plt.yticks(fontsize=14)    # Taille ticks y 

# Grille Lon/Lat
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()    
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter) 
# Define the yticks for latitude
ax.set_yticks(np.arange(-90,90,10), crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(axis='y', labelsize=14)
# Define the xticks for latitude
ax.set_xticks(np.arange(-180,180,20), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.tick_params(axis='x', labelsize=14)
ax.grid(linewidth=2, color='black', alpha=0.5, linestyle='') 

ax.set_title('SST bouées (obs) - {} '.format(date),fontsize=25)
figname = dir_fig + '{}_{}_{}_{}.png'.format(val2plot,instrument,exp,date)
print(figname)
plt.savefig(figname,dpi=100, format='png')
plt.clf()
plt.close("all")

# # Statistiques (data2plot)
# Moy_bouee=np.nanmean(data2plota)
# bouee_sst.append(Moy_bouee)
# Std_bouee=np.std(data2plota)
# bouee_std.append(Std_bouee)


# =============================================================================
# ## SST_bouees - Guess (ebauche)
# =============================================================================

print('*********** Ebauche')


data2plotb = odb_key_val_dic['col_8'][ind_flag]
lon       = odb_key_val_dic['col_1'][ind_flag]
lat       = odb_key_val_dic['col_2'][ind_flag]

val2plot='fg_depar_ass'
comp='bouees-ebauche'
instrument='arp'
exp_res="CMO"

   ## Colormap 
cmap_type = 'coolwarm' # le type de colormap
nb_of_bins = 9 # le nbre de couleurs dans la colormap
cmap = plt.get_cmap(cmap_type, nb_of_bins)
vmin, vmax = -2, 2
label_data="temp"

projection=ccrs.PlateCarree(central_longitude=-45)
figsize=(20,10)
fig,ax = plt.subplots(nrows=1,ncols=1,figsize = figsize,
                                    subplot_kw={'projection': projection},sharex=True,sharey=True)

im = ax.scatter(lon,lat,c=data2plotb, transform=ccrs.PlateCarree(),vmin = -2,vmax=2,cmap=cmap,s=50.0,marker="o")
ax.coastlines()
divider = make_axes_locatable(ax)
cax     = divider.append_axes("right", size="5%", pad=0.35, axes_class=plt.Axes)
bounds= [-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2]
cbar     = plt.colorbar(im,cax=cax,ticks=bounds,spacing='proportional',boundaries=[-2] + bounds + [2]) 
plt.ylabel('Biais des SST (K)', fontsize=20)
ticklabs = cbar.ax.get_yticklabels()
plt.yticks(fontsize=14)    # Taille ticks y 


# Grille Lon/lat
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()    
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
# Define the yticks for latitude
#ax.set_yticks(np.arange(-90,90,10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(20,71,10), crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(axis='y', labelsize=14)
# Define the xticks for latitude
ax.set_xticks(np.arange(-80,45,20), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.tick_params(axis='x', labelsize=14)
ax.grid(linewidth=2, color='black', alpha=0.5, linestyle='--') 

ax.set_title('SST bouées (obs) - SST ébauche {} ({}) - {}'.format(exp,exp_res,date),fontsize=25)
figname = dir_fig + '{}_{}_{}_{}.png'.format(comp,exp,val2plot,date)
print(figname)
plt.savefig(figname,dpi=100, format='png')


# Range data2plot
datamin   = min(data2plotb)
datamax   = max(data2plotb)
print('data2plot min est:',datamin)
print('data2plot max est:',datamax)
# # Statistiques (data2plot)
Moyenne=np.nanmean(data2plotb)
Ecart_type=np.nanstd(data2plotb)
Nb=len(data2plotb)
print('Moyenne bouées - ARP:',Moyenne)
print('Ecart-type bouées - ARP:',Ecart_type)
print('Nombre de points bouées - ARP:',Nb)
# # Tableau statistiques
# # fg
# T_stats_fg=np.array([[date,Moyenne,Ecart_type,Nb]])
# T_stats_date_fg=np.array([[date_2,Moyenne,Ecart_type,Nb]])
# T_data2plot_fg=np.array(data2plotb)
# print('Statistique à chaque réseau', T_stats_fg)
# # Liste et tableaux statistiques
# list_moy_fg.append(Moyenne)
# list_std_fg.append(Ecart_type)
# list_Nb_fg.append(Nb)
# list_date_fg.append(date_2)
# data2plot_fg=np.append(data2plot_fg,T_data2plot_fg,axis=0)
# T2_fg=np.append(T2_fg, T_stats_fg, axis=0)
# T_date2_fg=np.append(T_date2_fg, T_stats_date_fg, axis=0)
# #print(type(Stats_GDQV))
# #print(Stats_GDQV)



# =============================================================================
# ## SST bouees - analyse 
# =============================================================================



print('*********** Analyse')

data2plotc = odb_key_val_dic['col_9'][ind_flag]
lon       = odb_key_val_dic['col_1'][ind_flag]
lat       = odb_key_val_dic['col_2'][ind_flag]
 
val2plot='an_depar'
comp='bouees-analyse'
instrument='arp'
exp_res="CMO"

## Colormap 
cmap_type = 'coolwarm' # le type de colormap
nb_of_bins = 9 # le nbre de couleurs dans la colormap
cmap = plt.get_cmap(cmap_type, nb_of_bins)
vmin, vmax = -2, 2
label_data="temp"

projection=ccrs.PlateCarree(central_longitude=-45)
figsize=(20,10)
fig,ax = plt.subplots(nrows=1,ncols=1,figsize = figsize,
                            subplot_kw={'projection': projection},sharex=True,sharey=True)

im = ax.scatter(lon,lat,c=data2plotc, transform=ccrs.PlateCarree(),vmin =-2,vmax=2,cmap=cmap,s=50.0,marker="o")
ax.coastlines()
divider = make_axes_locatable(ax)
cax     = divider.append_axes("right", size="5%", pad=0.35, axes_class=plt.Axes)
bounds= [-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2]
cbar     = plt.colorbar(im,cax=cax,ticks=bounds,spacing='proportional',boundaries=[-2] + bounds + [2]) 
plt.ylabel('Biais des SST (K)', fontsize=20)
ticklabs = cbar.ax.get_yticklabels()
plt.yticks(fontsize=14)    # Taille ticks y 

# Plot des cartes BOUEES-ARP (GDQV)
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()    
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
# Define the yticks for latitude
#ax.set_yticks(np.arange(-90,90,10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(20,71,10), crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(axis='y', labelsize=14)
# Define the xticks for latitude
ax.set_xticks(np.arange(-80,45,20), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.tick_params(axis='x', labelsize=14)

ax.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')  
ax.set_title('SST bouées(obs) - SST analyse {} ({}) - {}'.format(exp,exp_res,date),fontsize=25)
figname = dir_fig + '{}_{}_{}_{}.png'.format(comp,exp,val2plot,date)
print(figname)
plt.savefig(figname,dpi=100, format='png')


# Range data2plot
datamin   = min(data2plotc)
datamax   = max(data2plotc)
print('data2plot min est:',datamin)
print('data2plot max est:',datamax)
# # Statistiques (data2plot)
Moyenne=np.nanmean(data2plotc)
Ecart_type=np.nanstd(data2plotc)
Nb=len(data2plotc)
print('Moyenne bouées - ARP:',Moyenne)
print('Ecart-type bouées - ARP:',Ecart_type)
print('Nombre de points bouées - ARP:',Nb)
#%%
## TEST WHOLE MAP
print('*********** WHOLE MAP')
# Importing the statistics module
import statistics
## SST bouees - analyse 
data2plotc = odb_key_val_dic['col_9']
lon       = odb_key_val_dic['col_1']
lat       = odb_key_val_dic['col_2']
 
val2plot='an_depar'
comp='bouees-analyse'
instrument='arp'
exp_res="CMO"

datamin   = min(data2plotc)
datamax   = max(data2plotc)
print('data2plot min est:',np.round(datamin,2))
print('data2plot max est:',np.round(datamax,2))
# # Statistiques (data2plot)
Moyenne=np.nanmean(data2plotc)
Ecart_type=np.nanstd(data2plotc)
Nb=len(data2plotc)
print('Moyenne bouées - ARP:',np.round(Moyenne,2))
print('Ecart-type bouées - ARP:',np.round(Ecart_type,2))
print('Nombre de points bouées - ARP:',Nb)

## Colormap 
cmap_type = 'coolwarm' # le type de colormap
nb_of_bins =11 # le nbre de couleurs dans la colormap
cmap = plt.get_cmap(cmap_type, nb_of_bins)
vmin, vmax = -2, 2
label_data="temp"

projection=ccrs.PlateCarree(central_longitude=-45)
figsize=(20,10)
fig,ax = plt.subplots(nrows=1,ncols=1,figsize = figsize,
                            subplot_kw={'projection': projection},sharex=True,sharey=True)

im = ax.scatter(lon,lat,c=data2plotc, transform=ccrs.PlateCarree(),vmin =-2,vmax=2,cmap=cmap,s=40.0,marker="o")
ax.coastlines()
divider = make_axes_locatable(ax)
cax     = divider.append_axes("right", size="5%", pad=0.35, axes_class=plt.Axes)
# bounds= [-5,-4.5,-4,-3.5,-3,-2.5,-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2, 2.5,3,3.5,4,4.5,5]
# bounds= [-5,-4,-3,-2,  -1, -0.5, -0.1, 0.1, 0.5, 1,2, 3,4,5]
bounds= [-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2]
cbar     = plt.colorbar(im,cax=cax,ticks=bounds,spacing='proportional',boundaries=[-2] + bounds + [2]) 
plt.ylabel('Biais (°C)', fontsize=20)
ticklabs = cbar.ax.get_yticklabels()
plt.yticks(fontsize=14)    # Taille ticks y 

# Plot des cartes BOUEES-ARP (GDQV)
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()    
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
# Define the yticks for latitude
#ax.set_yticks(np.arange(-90,90,10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(-90,90,10), crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(axis='y', labelsize=14)
# Define the xticks for latitude
ax.set_xticks(np.arange(-180,180,20), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.tick_params(axis='x', labelsize=14)

ax.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')  
ax.set_title('SST bouées(obs) - SST analyse {} ({}) - {}'.format(exp,exp_res,date),fontsize=25)
figname = dir_fig + '{}_{}_{}_{}.png'.format(comp,exp,val2plot,date)
print(figname)
plt.savefig(figname,dpi=100, format='png')

datamin   = min(data2plotc)
datamax   = max(data2plotc)
print('data2plot analyse min est:',datamin)
print('data2plot analyse max est:',datamax)

print()