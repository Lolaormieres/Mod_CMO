#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:19:23 2023

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






## PLOT CYCLE DIURNE - BOUEES
expe=['GN84','GN3C']
expe_name=['GN84']
SST_bouees=[]

### Ceate Dataframe
raw_data = {'date': [''],
            'lon': [''],
            'lat': [''],
        'expe': [''],
        'SST': ['']}
df_diurne_bouees= pd.DataFrame(raw_data, columns = ['date','lon','lat','expe','SST'])


## BOUES BOUCLE TEST
import datetime as dt
year=2022
month=9 # 4 pour Avril (et non 5 Mai)
day=1
hour=0
Delta=6
nbDate=17
Date_debut = dt.datetime(year, month, day, hour)

timedelta = dt.timedelta(0,0,0,0,0,Delta)

toto=[Date_debut+i*timedelta for i in range(nbDate)]

reseaux=""

# #Répertoire de travail
chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/GN84/'

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
            
    ind_flag = np.where((odb_key_val_dic['col_1']<5)&
    (odb_key_val_dic['col_1']>3)&
    (odb_key_val_dic['col_2']>37)&
    (odb_key_val_dic['col_2']<39))
        
    
    #data2plota = odb_key_val_dic['obsvalue'][ind_flag]
    data2plota =    odb_key_val_dic['col_7'][ind_flag]
    lon_bouees =    odb_key_val_dic['col_1'][ind_flag]
    lat_bouees =    odb_key_val_dic['col_2'][ind_flag]
    datamin   = 275
    datamax   = 300
    
    x_bouees_pt = data2plota[0] - 273.15
    print('len sst bouees', len(data2plota))
    # print('len lon bouees', len(lon_bouees))
    # print('len lat bouees', len(lat_bouees))
    
    print('sst bouees', x_bouees_pt)
    print('lon bouees', np.round(lon_bouees[0],2))
    print('lat bouees', np.round(lat_bouees[0],2))
    
    SST_bouees=np.append(SST_bouees,x_bouees_pt)
    print('SST append bouees', SST_bouees)

    # print(date_dt.strftime("%d%m%y%H%M"))
    date_3 = date_dt.strftime("%d%m%y%H%M")
    # date_3=int(date_3)
    # print(type(date_3))

    new_row = {'date':date_3, 'lon':lon_bouees[0], 'lat':lat_bouees[0], 'expe':expe[0],'SST': x_bouees_pt}
        
    # Use the loc method to add the new row to the DataFrame
    df_diurne_bouees.loc[len(df_diurne_bouees)] = new_row

# SUPRIMER 1er COL DF
df_diurne_bouees.drop(index=df_diurne_bouees.index[0], axis=0, inplace=True) ##

#%%

lon_atm= 4.1
lat_atm= 38.1

# lon bouees 3.16
# lat bouees 39.33

# CARTE BOUEEDS - ILES BALEARES

 
## Colormap 
cmap_type = 'Blues' # le type de colormap
nb_of_bins = 9 # le nbre de couleurs dans la colormap
cmap = plt.get_cmap(cmap_type, nb_of_bins)
vmin, vmax = -2, 2
label_data="temp"

projection=ccrs.PlateCarree(central_longitude=-45)
figsize=(20,10)
fig,ax = plt.subplots(nrows=1,ncols=1,figsize = figsize,
                        subplot_kw={'projection': projection},sharex=True,sharey=True)
# im = ax.scatter(lon_atm,lat_atm,c=data2plota[0]+6, transform=ccrs.PlateCarree(),vmin =-2,vmax=2,cmap=cmap,s=70.0,marker="o")
im = ax.scatter(lon_bouees[0],lat_bouees[0],c=data2plota[0], transform=ccrs.PlateCarree(),vmin =-2,vmax=2,cmap=cmap,s=150.0,marker="o")
#im = ax.scatter(lon_atm,lat_atm,c=data2plota[0]+6, transform=ccrs.PlateCarree(),vmin =-2,vmax=2,cmap=cmap,s=70.0,marker="o")
ax.coastlines()
divider = make_axes_locatable(ax)
cax     = divider.append_axes("right", size="5%", pad=0.5, axes_class=plt.Axes)
bounds= [-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2]
cbar     = plt.colorbar(im,cax=cax,ticks=bounds,spacing='proportional',boundaries=[-2] + bounds + [2]) 
# cbar_ax = fig.add_axes()
   
# cax_ax.tick_params(labelsize=18)

# cbar.ax.set_title('Bias (°C)',fontsize=19)
# cbar.ax.get_yticklabels(fontsize=16)
cbar.ax.tick_params(labelsize=26)

plt.ylabel('Bias (°C)', fontsize=20)
ticklabs = cbar.ax.get_yticklabels()
plt.yticks(fontsize=22)    # Taille ticks y 

# Plot des cartes BOUEES-ARP (GDQV)
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()    
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
# Define the yticks for latitude
#ax.set_yticks(np.arange(-90,90,10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(20,71,10), crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(axis='y', labelsize=18)
# Define the xticks for latitude
ax.set_xticks(np.arange(-35,46,10), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.tick_params(axis='x', labelsize=18)

ax.grid(linewidth=1, color='black', alpha=0.5, linestyle='--')  
# ax.set_title('SST bouées(obs) - SST analyse {} ({}) - {}'.format(exp,exp_res,date),fontsize=25)
ax.set_title('Analised SST ({}) - SST Buoys, {}'.format(expe_name,date_dt),fontsize=25)
plt.show()
# figname = dir_fig + '{}_{}_{}_{}.png'.format(zone,date,expe,comp)
# plt.savefig(figname,dpi=100, format='png')