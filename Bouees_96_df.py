#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 17:13:52 2023

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
pathres="/home/ormieresl/Routines/bouees/DF_bouees_96h"

# experience to evaluate
expe=['GN84','GO4A']
# reference
expe_mer=['GO48']
compar='prev-mer'+expe_mer[0]

# Levels numbers
nb_level=12

# Start date
date_init = datetime(2022, 9,1 , 0)
# Days numbers
nb=2
nbDate=nb # buyos dates

# Term
ech_fin=102
ech_ini=0
ech_step=6

# Variables names
name_sst=['SFX.SST']
name=['SFX.TEMP_OC'+str(i) for i in range(1,nb_level+1)]
name_sal=['SFX.SALT_OC'+str(i) for i in range(1,nb_level+1)]
name_mer=['SURF.THETAO'+str(i) for i in range(1,nb_level+1)]
name_sal_mer=['SURF.SALINO'+str(i) for i in range(1,nb_level+1)]
name_u=['SFX.UCUR_OC'+str(i) for i in range(1,nb_level+1)]
name_v=['SFX.VCUR_OC'+str(i) for i in range(1,nb_level+1)]

## Level file
filename = 'level2.txt'
levelfic = np.loadtxt(filename, delimiter=',', dtype=float)
# print('level', levelfic)
# print(type(levelfic))
# print(np.shape(levelfic))
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

# =============================================================================
# RECUPERER OBS VALUE GO4A
# =============================================================================

## Lecture fichier Netcdf
# expe=['GO4A']
expe_name='Xp.assim.buoys'
comp='analyse-bouees'
dir_fig='Cartes_bouees_v1/'
zone='eurat'

biais_obs=[[]]
chemin='/home/ormieresl/Scripts/'
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
cterm_for='96'
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
########################################################
# =============================================================================
# ## Lectur ficher fa
# =============================================================================
#LECTUR FILES
def get_file(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,
          'local':'tmpbouees.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpbouees.fic',openmode='r')
    os.remove('tmpbouees.fic')
    return r
#%%
# =============================================================================
# ## Time series: P96 - buyos --> CSV
# =============================================================================
### Ceate Dataframe FINAL
raw_data =      {'date_valid': [''],
                 'date_xp': [''],
                 'expe': [''],
                 'level': [''],
                 'ech': [''],
                 'biais_obs': [''],
                 'eqm': [''],
                 'Nb_bouees': [''],
                 'diff_Nb':['']}
df_bouees_for = pd.DataFrame(raw_data, columns = ['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees','diff_Nb'])


# Date
# -----------------------------------------------------------------------------

# dt = timedelta(days=1)
# date_reseau = date_init - dt
# date_valid=date_init-timedelta(hours=ech_step)
# period = timedelta(days=nb)
# date_end = date_reseau+period


for e in range(len(expe)):
    print(expe[e]) 
    import datetime as dt
    date_xp = date_init
    deltat = timedelta(days=1)
    period = timedelta(days=nb)
    date_end = date_init+period
    for t in range(nb):
        print('date_xp',date_xp)
        biais_obs=[]
        biais_obs=np.array(biais_obs)
        ## DF daily difference
        raw_data={'biais_obs':[''],'lon':[''],'lat':['']}
        df_diff_day = pd.DataFrame(raw_data, columns = ['biais_obs','lon','lat'])
        
        # =============================================================================
        # Buyos Loop
        # =============================================================================
        year=2022
        month=9 
        day=5
        hour=0
        Delta=24
        #nbDate=1
        Date_debut = dt.datetime(year, month, day, hour)
        timedelta = dt.timedelta(0,0,0,0,0,Delta)
        toto=[Date_debut+i*timedelta for i in range(nb)]
        reseaux=""
        ##XP BOUCLE TEST
        year=2022
        month=9
        day=1
        hour=0
        Delta=24
        #nbDate=1
        Date_debut_xp = dt.datetime(year, month, day, hour)
        timedelta = dt.timedelta(0,0,0,0,0,Delta)
        toto_xp=[Date_debut_xp+i*timedelta for i in range(nb)]
        reseaux=""
        
        #Data repertory EXPE
        chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(expe[e])+'/'
        
        #BOUCLE read files
        # cterm_for= ['0','6','12','18','0','6','12','18','0','6','12','18','0','6','12','18','0']
        # for date_dt in toto:
        date_dt=toto[t]
        date_dt_expe=toto_xp[t]
        
        print('date_valid',date_dt)
        print('date_xp',date_dt_expe)
        
        #Date - Repertory
        date = '%4.4i%2.2i%2.2iT%2.2i00A'%(date_dt.year,\
        date_dt.month,\
        date_dt.day,\
        date_dt.hour)
        print(date)
        date_2 = '%2.2i/%2.2i/%2.2i'%(date_dt.year,\
        date_dt.month,\
        date_dt.day)
        print(date_2)
        date_expe = '%2.2i/%2.2i/%2.2i'%(date_dt_expe.year,\
        date_dt_expe.month,\
        date_dt_expe.day) # date_dt.day-4,\ works but issue when we do over a month period
        date_valid=date_2
        date_xp = date_xp + deltat
    
        print('clean date  valid***********',date_valid)
        print('clean date  xp**********',date_xp)
         
        fic_nc = chemin_nc +'{}/surfan/odb-ecma.canari.surf.nc'.format(date)
        # print(fic_nc)
        if not os.path.isfile(fic_nc):
            # print('pas de fichier '+fic_nc)
            print('pas de fichier ')
        #print(len(fic_nc))   
    
        if os.path.isfile(fic_nc):
        #Ouverture du fichier netcdf au format dataset
            # print('Ouverture du fichier',fic_nc)
            print('Ouverture du fichier')
            myfile = Dataset(fic_nc,'r')
        #print(myfile.variables) #structure fichier
    
        #Pour donner des noms explicites aux différentes variables :
        odb_key_val_dic = {}
        for variables in myfile.variables.keys():
        #instrument_int = myfile[variables].odb_name[:myfile[variables].odb_name.index('@')]
            odb_key_val_dic[variables] = myfile[variables][:]
                #print(variables,'->',instrument_int)
                    
        #Flag zone
        ind_flag = np.where((odb_key_val_dic['col_1']<45)&
        (odb_key_val_dic['col_1']>-35)&
        (odb_key_val_dic['col_2']>20)&
        (odb_key_val_dic['col_2']<72))    
        
        #an-depar
        data2plotc = -odb_key_val_dic['col_9'][ind_flag] ## 9 an-depar (bouees-analyse), 8 fg-depar
        lon       = odb_key_val_dic['col_1'][ind_flag]
        lat       = odb_key_val_dic['col_2'][ind_flag]
        id_bouee  = odb_key_val_dic['col_5'][ind_flag]
        
        #obs value
        data2plota = odb_key_val_dic['col_7'][ind_flag] ## 7 obsvalue
        lon_2      = odb_key_val_dic['col_1'][ind_flag]
        lat_2      = odb_key_val_dic['col_2'][ind_flag]
        id_bouee_obs  = odb_key_val_dic['col_5'][ind_flag]
        date      = odb_key_val_dic['col_3'][ind_flag]
        time      = odb_key_val_dic['col_4'][ind_flag]
        
        #DF buyos data from nc files
        a=pd.DataFrame(lon_2)
        b=pd.DataFrame(lat_2)
        c=pd.DataFrame(data2plota-273.15)
        #ID=pd.DataFrame(id_bouee_obs)
        date=pd.DataFrame(time)
        temps=pd.DataFrame(time)
        
        #DATAFRAME BUYOS -- meme bouées pas même heure
        df_bouees=pd.concat([a,b,c,date,temps],axis=1)
        df_bouees.columns =['lon', 'lat', 'obs','date','time']
        #obs_sst=data2plotc+data2plota
        obs_sst=data2plota
        #dates = pd.date_range('1/1/2000', periods=8)
        
        #Buyos numbers
        Nb_obs_assim=len(data2plotc)
        Nb_obs=len(data2plota)
        Nb_diff=Nb_obs-Nb_obs_assim
        
        #test lon at one grid point
        lon_test=lon[0]
        lat_test=lat[0]
        # print('lon test',lon_test)
        # print('lat test', lat_test)
       

        # =============================================================================
        # XP
        # =============================================================================

        #Read Xp data
        r_for=get_file(expe[e],date_xp,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
    
        #Loop on buoys lon/lat, line after line of Buyos DF
        for i in range(Nb_obs):
            #for i in lon:
            #    for j in lat:
                    # print('i***********************************************',i)
                    # print('j***********************************************',j)
            # =============================================================================
            # XP
            # =============================================================================
            
            X=r_for.geometry.get_lonlat_grid()
            lon_expe=X[0]
            lat_expe=X[1]
            
            x=r_for.readfield(name[0])
                        
            SST=x.extract_point(df_bouees.lon[i],df_bouees.lat[i]) 
            # print(SST) 
            # print('lon',lon_test)
            # print('lat',lat_test)       
            
            # =============================================================================
            #Calcul diff --> XP - OBS
            # =============================================================================
            # for it in range(len(obs_sst)):
            biais=SST-df_bouees.obs[i]
            biais=np.array(biais.data)
            new_row = {'biais_obs':biais,'lon':df_bouees.lon[i], 'lat':df_bouees.lat[i]}    
            # print(new_row)
            df_diff_day.loc[len(df_diff_day)] = new_row 
        df_diff_day.drop(index=df_diff_day.index[0], axis=0, inplace=True)
        # print('df_obs_diff, avant mean',df_diff_day)
         
            
        # biais_obs=np.append(biais_obs,biais.data)
        # biais_obs_mean=np.mean(biais_obs) ## Marche pas
        biais_obs_mean=df_diff_day['biais_obs'].mean()
        biais_obs_std=df_diff_day['biais_obs'].std()
        # print('Daily Mean biais: ',biais_obs_mean)
                
                    # np.hstack # add horizontal
                    # np.vstack # add vertical
            
        # date_expe = date_expe.strftime("%y/%m/%d")
        new_row = {'date_valid':date_valid,'date_xp':date_expe, 'expe':expe[e], 'var':name[0],'level':level[0], 'ech':cterm_for,'biais_obs': biais_obs_mean,'eqm': biais_obs_std,'Nb_bouees':Nb_obs ,'diff_Nb':Nb_diff}    
        # print(new_row)
        df_bouees_for.loc[len(df_bouees_for)] = new_row  # Use the loc method to add the new row to the DataFrame
print('df_bouees_for final:', df_bouees_for)
df_bouees_for.drop(index=df_bouees_for.index[0], axis=0, inplace=True) # SUPRIMER 1er COL DF

 
df_bouees_for.to_csv(pathres+"/"+str(zone)+'_P96-bouees_septembre.csv')



# df_bouees_for['biais_obs'].mean()
#%%
# =============================================================================
#Plot time series of SST difference
# =============================================================================
fig, ax = plt.subplots(nrows = 1,ncols = 1, sharex = False, figsize=(8, 5),dpi=100) 
N=2
L=len(df_bouees_for.date_valid)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_bouees_for.date_valid[0:len(df_bouees_for.date_valid):N]   

ax=df_bouees_for.plot(x='date_valid', y='biais_obs', kind='scatter',color='orange', ax=ax) ## kind: scatter, line
ax.set_xticks([])
ax.set_ylabel('Biais (°C)',fontsize=16)
xlabels = df_bouees_for.date_valid[0:len(df_bouees_for.date_valid):N]
ax.set_xticks(xind,labels=xlabels,fontsize=13,rotation=45)
ax.set_xlabel('Date of obs',fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)
ax.set_title('Difference of SST between CMO mod. and buyos',fontsize=20)
    

#%%
# =============================================================================
#Plot Map of bbiais (put in loop, to have one map per day)
# =============================================================================
# ## Plot Data2plotc - an-depar

# ## Colormap 
# cmap_type = 'coolwarm' # le type de colormap
# nb_of_bins = 9 # le nbre de couleurs dans la colormap
# cmap = plt.get_cmap(cmap_type, nb_of_bins)
# vmin, vmax = -2, 2
# label_data="temp"

# projection=ccrs.PlateCarree(central_longitude=-45)
# figsize=(20,10)
# fig,ax = plt.subplots(nrows=1,ncols=1,figsize = figsize,
#                         subplot_kw={'projection': projection},sharex=True,sharey=True)

# im = ax.scatter(lon,lat,c=-data2plotc, transform=ccrs.PlateCarree(),vmin =-2,vmax=2,cmap=cmap,s=50.0,marker="o")
# ax.coastlines()
# divider = make_axes_locatable(ax)
# cax     = divider.append_axes("right", size="5%", pad=0.35, axes_class=plt.Axes)
# bounds= [-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2]
# cbar     = plt.colorbar(im,cax=cax,ticks=bounds,spacing='proportional',boundaries=[-2] + bounds + [2]) 
# plt.ylabel('Biais des SST (K)', fontsize=20)
# ticklabs = cbar.ax.get_yticklabels()
# plt.yticks(fontsize=14)    # Taille ticks y 

# # Plot des cartes BOUEES-ARP (GDQV)
# lon_formatter = cticker.LongitudeFormatter()
# lat_formatter = cticker.LatitudeFormatter()    
# ax.xaxis.set_major_formatter(lon_formatter)
# ax.yaxis.set_major_formatter(lat_formatter)
# # Define the yticks for latitude
# ax.set_yticks(np.arange(20,71,10), crs=ccrs.PlateCarree())
# ax.yaxis.set_major_formatter(lat_formatter)
# ax.tick_params(axis='y', labelsize=14)
# # Define the xticks for latitude
# ax.set_xticks(np.arange(-35,46,10), crs=ccrs.PlateCarree())
# ax.xaxis.set_major_formatter(lon_formatter)
# ax.tick_params(axis='x', labelsize=14)

# ax.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')  
# ax.set_title('SST analyse  {} - SST bouées(obs)  - {}'.format(expe,date),fontsize=25)
# # figname = dir_fig + '{}_{}_{}_{}.png'.format(zone,date,expe,comp)
# # plt.savefig(figname,dpi=100, format='png')

 
# ## Plot Obs value

# ## Colormap 
# cmap_type = 'plasma' # le type de colormap
# nb_of_bins = 9 # le nbre de couleurs dans la colormap
# cmap = plt.get_cmap(cmap_type, nb_of_bins)
# vmin, vmax = 270, 300
# label_data="temp"

# projection=ccrs.PlateCarree(central_longitude=-45)
# figsize=(20,10)
# fig,ax = plt.subplots(nrows=1,ncols=1,figsize = figsize,
#                         subplot_kw={'projection': projection},sharex=True,sharey=True)

# im = ax.scatter(lon,lat,c=obs_sst, transform=ccrs.PlateCarree(),vmin =270,vmax=300,cmap=cmap,s=50.0,marker="o")
# ax.coastlines()
# divider = make_axes_locatable(ax)
# cax     = divider.append_axes("right", size="5%", pad=0.35, axes_class=plt.Axes)
# #np.arange(200,310,10)
# # bounds= [200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300]
# bounds= [ 270, 280, 290, 300]
# cbar     = plt.colorbar(im,cax=cax,ticks=bounds,spacing='proportional',boundaries=[270] + bounds + [300]) 
# plt.ylabel('Biais des SST (K)', fontsize=20)
# ticklabs = cbar.ax.get_yticklabels()
# plt.yticks(fontsize=14)    # Taille ticks y 

# lon_formatter = cticker.LongitudeFormatter()
# lat_formatter = cticker.LatitudeFormatter()    
# ax.xaxis.set_major_formatter(lon_formatter)
# ax.yaxis.set_major_formatter(lat_formatter)

# # Define the yticks for latitude
# ax.set_yticks(np.arange(20,71,10), crs=ccrs.PlateCarree())
# ax.yaxis.set_major_formatter(lat_formatter)
# ax.tick_params(axis='y', labelsize=14)
# # Define the xticks for latitude
# ax.set_xticks(np.arange(-35,46,10), crs=ccrs.PlateCarree())
# ax.xaxis.set_major_formatter(lon_formatter)
# ax.tick_params(axis='x', labelsize=14)
# ax.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')  

# ax.set_title('SST analyse  {} - SST bouées(obs)  - {}'.format(expe,date),fontsize=25)
# # figname = dir_fig + '{}_{}_{}_{}.png'.format(zone,date,expe,comp)
# # plt.savefig(figname,dpi=100, format='png')

