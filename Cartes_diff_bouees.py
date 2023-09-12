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
expe='GO4A'
comp='analyse-bouees'
dir_fig='Cartes_bouees/'
zone='eurat'
## Date - exp
nb=1
date_init = datetime(2022, 10,20, 0)
date = date_init
dt = timedelta(days=1)
date_valid = date+dt
period = timedelta(days=nb)
date_end = date+period

## PLOT CYCLE DIURNE - BOUEES

SST_bouees=[]

### Ceate Dataframe
raw_data = {'date': [''],
            'lon': [''],
            'lat': [''],
        'expe': [''],
        'biais': [''],
        'id_bouee':['']}
df_bouees= pd.DataFrame(raw_data, columns = ['date','lon','lat','expe','biais','id_bouee'])

## BOUES BOUCLE TEST
import datetime as dt
year=2022
month=10 # 4 pour Avril (et non 5 Mai)
day=7
hour=0
Delta=0
nbDate=1
Date_debut = dt.datetime(year, month, day, hour)

timedelta = dt.timedelta(0,0,0,0,0,Delta)

toto=[Date_debut+i*timedelta for i in range(nbDate)]

reseaux=""

# #Répertoire de travail
chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(expe)+'/'

# BOUCLE lecture des fichiers
# cterm_for= ['0','6','12','18','0','6','12','18','0','6','12','18','0','6','12','18','0']
for date_dt in toto:
    print(date_dt)
    date = '%4.4i%2.2i%2.2iT%2.2i00A'%(date_dt.year,\
    date_dt.month,\
    date_dt.day,\
    date_dt.hour)
    print(date)
    # print(toto)
    date_2 = '%2.2i/%2.2i/%2.2i %2.2i'%(date_dt.year,\
    date_dt.month,\
    date_dt.day,\
    date_dt.hour)
    print(date_2)
    date_3 = '%2.2i%2.2i%2.2i%2.2i00'%(date_dt.year,\
    date_dt.month,\
    date_dt.day,\
    date_dt.hour)
    print(date_3)
    

    fic_nc = chemin_nc +'{}/surfan/odb-ecma.canari.surf.nc'.format(date)
    print(fic_nc)
    if not os.path.isfile(fic_nc):
        print('pas de fichier '+fic_nc)
    #print(len(fic_nc))   

    if os.path.isfile(fic_nc):
    #Ouverture du fichier netcdf au format dataset
        print('Ouverture du fichier',fic_nc)
        myfile = Dataset(fic_nc,'r')
    #print(myfile.variables) #structure fichier

    #Pour donner des noms explicites aux différentes variables :
    odb_key_val_dic = {}
    for variables in myfile.variables.keys():
    #instrument_int = myfile[variables].odb_name[:myfile[variables].odb_name.index('@')]
        odb_key_val_dic[variables] = myfile[variables][:]
            #print(variables,'->',instrument_int)
                
            
    # ind_flag = ((odb_key_val_dic['col_5'].data==[b' ', b'6', b'2', b'0', b'4', b'5', b'9', b'4']))
    
    ind_flag = np.where((odb_key_val_dic['col_1']<-130)&
    (odb_key_val_dic['col_1']<60)&
    (odb_key_val_dic['col_2']>50)&
    (odb_key_val_dic['col_2']<72)&
    (odb_key_val_dic['col_10'] ==1))    
    
    # ind_flag = np.where((odb_key_val_dic['col_1']<45)&
    # (odb_key_val_dic['col_1']>-35)&
    # (odb_key_val_dic['col_2']>20)&
    # (odb_key_val_dic['col_2']<72)&
    # (odb_key_val_dic['col_10'] ==1))    
    
    # ind_flag = np.where((odb_key_val_dic['col_1']<45)&
    # (odb_key_val_dic['col_1']>-35)&
    # (odb_key_val_dic['col_2']>20)&
    # (odb_key_val_dic['col_2']<72)&
    # (odb_key_val_dic['col_5'].data==[b' ', b'6', b'2', b'0', b'4', b'5', b'9', b'4']))
    
    
    data2plotc = -odb_key_val_dic['col_8'][ind_flag]
    lon       = odb_key_val_dic['col_1'][ind_flag]
    lat       = odb_key_val_dic['col_2'][ind_flag]
    id_bouee  = odb_key_val_dic['col_5'][ind_flag]
    
#     #Creating an empty list to to append the SICs for each netCDF file
#     b = []
# #Looping through the entire files and creating a variable for each SIC
#     for i in range(len(nbDate)):
#         exec (f'df{i} = Dataset(w[i])')
#         exec (f"sea_ice_concentration{i} = df{i}.variables['z'][:]")
#         exec (f"sic{i} = np.flip(sea_ice_concentration{i},0)")
#         exec(f'b.append(sic{i})')
   
        
   
    

    
    print('Biais moy', np.mean(data2plotc))
    print('min biais', min(data2plotc))
    print('max biais', max(data2plotc))
    print('Nb bouées', len(data2plotc))
    
    new_row = {'date':date, 'lon':lon, 'lat':lat, 'expe':expe,'biais':data2plotc,'id_bouee':id_bouee}
        
    df_bouees.loc[len(df_bouees)] = new_row
    
    # # ## Sauvegarde dataframe 
    # df_bouees.to_csv(dir_fig+'/analyse_bouees'+str(expe)+'.csv')   
    
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
    ax.set_yticks(np.arange(50,71,10), crs=ccrs.PlateCarree())
    ax.yaxis.set_major_formatter(lat_formatter)
    # ax.tick_params(axis='y', labelsize=14)
    # Define the xticks for latitude
    # ax.set_xticks(np.arange(-200,46,10), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.tick_params(axis='x', labelsize=14)

    ax.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')  
    # ax.set_title('SST bouées(obs) - SST analyse {} ({}) - {}'.format(exp,exp_res,date),fontsize=25)
    ax.set_title('SST analyse  {} - SST bouées(obs)  - {}'.format(expe,date),fontsize=25)
    figname = dir_fig + '{}_{}_{}_{}.png'.format(zone,date,expe,comp)
    # plt.savefig(figname,dpi=100, format='png')
    
    