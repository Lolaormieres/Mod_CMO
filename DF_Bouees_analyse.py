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
dir_fig='/cnrm/recyf/Data/users/ormieresl/plot_evol_sst_bouees/'
# experience to evaluate
expe=['GN84','GO4A']
# expe=['GKCJ','GO48']


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

cterm_for='96'
zone =  'nordat'

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
cterm_for='6'
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
          'local':'tmpbouees_test.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpbouees_test.fic',openmode='r')
    os.remove('tmpbouees_test.fic')
    return r

#%%
# =============================================================================
# ## Time series: P96 - buyos --> CSV
# =============================================================================
#DF / period   --> Initialised DataFrame final

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
import datetime as dt
# Buyos Loop
year=2022
month=9 
day=2
hour=0
Delta=24
#nbDate=1
Date_debut = dt.datetime(year, month, day, hour)
Date_debut_valid = dt.datetime(year, month, day, hour) # Date_valid = Date XP
timedelta = dt.timedelta(0,0,0,0,0,Delta)
toto=[Date_debut+i*timedelta for i in range(nbDate)]
toto_valid=[Date_debut_valid+i*timedelta for i in range(nbDate)]

for e in range(len(expe)):
    print(expe[e])

    date_xp = date_init - deltat # Reinitialise Date between each xp
    date_end = date_xp+period
    
    #Data repertory EXPE
    chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(expe[e])+'/'
    
    for t in range(nb):
        date_xp = date_xp + deltat
          # Loop on expe files timestep = 1 day
        # print('date_xp',date_xp)
        
        #DF / Daily   --> Initialised DataFrame in loop
        raw_data={'biais_obs':[''],'lon':[''],'lat':['']}
        df_diff_day = pd.DataFrame(raw_data, columns = ['biais_obs','lon','lat'])
        
        date_dt=toto[t]
        date_dt_valid=toto_valid[t]
        # print('date_xp',date_dt)
        # print('date_valid',date_dt_valid)
        
        #Date - Repertory - Buyos
        date_buyos_dir = '%4.4i%2.2i%2.2iT%2.2i00A'%(date_dt_valid.year,\
        date_dt_valid.month,\
        date_dt_valid.day,\
        date_dt_valid.hour)
        print('Date fic NC',date_buyos_dir)
        
        #Date Valid Clean
        date_valid = '%2.2i/%2.2i/%2.2i'%(date_dt_valid.year,\
        date_dt_valid.month,\
        date_dt_valid.day)
        
        #Date Clean
        date_expe = '%2.2i/%2.2i/%2.2i'%(date_dt.year,\
        date_dt.month,\
        date_dt.day) 
         
        # DATE for DataFrame
        # print('clean date VALID',date_valid)
        # print('clean date XP',date_expe)
          
        #Read NetCDF
        fic_nc = chemin_nc +'{}/surfan/odb-ecma.canari.surf.nc'.format(date_buyos_dir)
        
        #Test on files
        if not os.path.isfile(fic_nc):
            print('pas de fichier '+fic_nc)
            miss_file=len(fic_nc)           # Variables with number of missing NC files 
        if os.path.isfile(fic_nc):
            print('Ouverture du fichier')   # print('Ouverture du fichier',fic_nc)
            myfile = Dataset(fic_nc,'r')    #Ouverture du fichier netcdf au format dataset     #print(myfile.variables) #structure fichier
        
        #Give explicites name at different variables of Nc files:
        odb_key_val_dic = {}
        for variables in myfile.variables.keys():
            odb_key_val_dic[variables] = myfile[variables][:]
    
        
        
        # Flag zone --> Zone (Do dict) N.atl
        ind_flag = np.where((odb_key_val_dic['col_1'] < 0) &
                            (odb_key_val_dic['col_1'] > -80) &
                            (odb_key_val_dic['col_2'] > 20) &
                            (odb_key_val_dic['col_2'] < 70) &
                            (odb_key_val_dic['col_10'] ==1))
        

        # Flag zone --> Zone (Do dict) eurat
        ind_flag_off = np.where((odb_key_val_dic['col_1'] < 45) &
                                (odb_key_val_dic['col_1'] > -35) &
                                (odb_key_val_dic['col_2'] > 20) &
                                (odb_key_val_dic['col_2'] < 72) &     
                                (odb_key_val_dic['col_10'] ==1))

        # Flag zone --> Zone (Do dict) tropiques
        ind_flag_off = np.where((odb_key_val_dic['col_2'] < 20) &
                                (odb_key_val_dic['col_2'] > -20) &
                                (odb_key_val_dic['col_10'] ==1))

        # Flag zone --> Zone (Do dict) hm n
        ind_flag_off = np.where((odb_key_val_dic['col_2'] > 20) & 
                                (odb_key_val_dic['col_10'] ==1))

        # Flag zone --> Zone (Do dict) hm n
        ind_flag_off = np.where((odb_key_val_dic['col_2'] < -20)
                                & (odb_key_val_dic['col_10'] ==1))

        #an-depar
        data2plotc = -odb_key_val_dic['col_9'][ind_flag]    ## 9 an-depar (bouees-analyse), 8 fg-depar
        lon       = odb_key_val_dic['col_1'][ind_flag]
        lat       = odb_key_val_dic['col_2'][ind_flag]
        id_bouee  = odb_key_val_dic['col_5'][ind_flag]
        
        #obs value
        data2plota = odb_key_val_dic['col_7'][ind_flag]     ## 7 obsvalue
        lon_2      = odb_key_val_dic['col_1'][ind_flag]
        lat_2      = odb_key_val_dic['col_2'][ind_flag]
        id_obs  = odb_key_val_dic['col_5'][ind_flag]
        date      = odb_key_val_dic['col_3'][ind_flag]      # in secondes
        time      = odb_key_val_dic['col_4'][ind_flag]      # in secondes
        
        #DF buyos data from nc files
        a=pd.DataFrame(lon)
        b=pd.DataFrame(lat)
        c=pd.DataFrame(data2plota-273.15)
        date=pd.DataFrame(time)                             # secondes
        temps=pd.DataFrame(time)                            # secondes
        
        
        #DATAFRAME BUYOS -- meme bouées pas même heure
        df_bouees=pd.concat([a,b,c,date,temps],axis=1)
        df_bouees.columns =['lon', 'lat', 'obs','date','time']
        obs_sst=data2plota

        #Buyos numbers
        Nb_obs_assim=len(data2plotc)
        Nb_obs=len(data2plota)
        Nb_diff=Nb_obs-Nb_obs_assim  #Diff between nubers of obs tot and numbers of obs assim
        
        #Read Xp data
        # =============================================================================
        #r_for=get_file(expe[e],date_xp,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
        r=get_file(expe[e],date_xp,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
        X=r.geometry.get_lonlat_grid()
        lon_expe=X[0]
        lat_expe=X[1]
        x=r.readfield(name[0])
        
        #Loop on dataframe of obs
        # =============================================================================
        for i in range(Nb_obs):
            SST=x.extract_point(df_bouees.lon[i],df_bouees.lat[i])   #Extract data on lon/lat
            # print(SST)
            biais=SST.data-df_bouees.obs[i]
            
            new_row1 = {'biais_obs':biais,'lon':df_bouees.lon[i], 'lat':df_bouees.lat[i]} 
            df_diff_day.loc[len(df_diff_day)] = new_row1 
            
        df_diff_day.drop(index=df_diff_day.index[0], axis=0, inplace=True) # delete firt row
        df_diff_day.drop(df_diff_day.loc[df_diff_day['biais_obs']>1.5e+100].index, inplace=True) # delete Nan values     
         
        # Compute mean diff and std over wone, eache day
        biais_obs_mean=df_diff_day['biais_obs'].mean()
        biais_obs_std=df_diff_day['biais_obs'].std()
        
        new_row = {'date_valid':date_valid,'date_xp':date_expe, 'expe':expe[e], 'var':name[0],'level':level[0], 'ech':cterm_for,'biais_obs': biais_obs_mean,'eqm': biais_obs_std,'Nb_bouees':Nb_obs ,'diff_Nb':Nb_diff}    
        df_bouees_for.loc[len(df_bouees_for)] = new_row 
    
#Save DataFrame
df_bouees_for.drop(index=df_bouees_for.index[0], axis=0, inplace=True) # delete first row

df_bouees_for.to_csv(pathres+str(zone)+'_A'+str(hour)+'h-bouees_flag.csv') 


