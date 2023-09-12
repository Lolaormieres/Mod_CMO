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
#pathres='/home/ormieresl/Routines/DF_bouees_ech_vf/DF_bouees_ech_vf_maskcorr/'
pathres = '/home/ormieresl/Routines/DF_BOUEES_mask/longerperiod/'
dir_fig='/cnrm/recyf/Data/users/ormieresl/plot_ech_bouees/'


# experience to evaluate
expe=['GN84','GO48','GKCJ']
#expe=['GOJQ']
# expe=['GKCJ','GO48']


# reference
expe_mer=['GO48']
compar='prev-mer'+expe_mer[0]

# Levels numbers
nb_level=1

# Start date
date_init = datetime(2022,8,2, 0)
deltat = timedelta(days=1)
ech = [hour for hour in np.arange(0, 102, 6)]

# Days numbers
nb=180
nbDate=nb # buyos dates
period = timedelta(days=nb)

# échéances
ech_fin=102
ech_ini=0
ech_step=6

#zone =  'nordat'

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
color_GL55='#648FFF'  
color_GKOR='#785EF0'  
color_GKPH='red'  
color_GM6E='#FE6100'  
color_GMOT='#009E73'
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
fill_for='surf'

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
          'local':'tmpb_ech_mask.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpb_ech_mask.fic',openmode='r')
    os.remove('tmpb_ech_mask.fic')
    return r

#%%
# =============================================================================
# ## Name variables & dates
# =============================================================================
#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [-180,180,-20,20],"hn20":[-180,180,20,90],"hs20":[-180,180,-90,-20],"glob":[-180,180,-90,90],"med": [-3,16,30,44]}

# =====================date_valid.hour=========================================
# ### Previ - Mercator 
# =============================================================================
## Prev - mer 
#DF / period   --> Initialised DataFrame finaldate_reseau = date_init - dt
raw_data =      {'date_valid': [''],
                 'date_xp': [''],
                 'expe': [''],
                 'level': [''],
                 'ech': [''],
                 'biais_obs': [''],
                 'eqm': [''],
                 'Nb_bouees_tot': [''],
                 'Nb_bouees_an':[''],
                 'diff_Nb':['']}

####################
# start
####################
for e in range(len(expe)):
    df_bouees_final={}
    for zon in zones.keys():
        df_bouees_final[zon] = pd.DataFrame(raw_data, columns = ['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees_tot','Nb_bouees_an','diff_Nb'])
    
    chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(expe[e])+'/'
    dt = timedelta(days=1)
    date_reseau = date_init - dt

    date_valid=date_init-timedelta(hours=ech_step)

    period = timedelta(days=nb)
    date_end = date_reseau+period
    
             
    #Mask for buyos data zone
    for key in zones.keys():
        print('key',key)
        print('date_reseauuuuuuuuuuuuuuuuuuu',date_reseau)
        date_reseau = date_init - dt            #### reinitialise date ici!!!
        for t in range(nb):
            date_reseau = date_reseau+dt
            cterm_for= range(ech_ini,ech_fin,ech_step)
            
            for h in range(len(cterm_for)):
               print('ech',h)
               dt2=timedelta(hours=cterm_for[h])
               date_valid = date_reseau+dt2
                
               #DF / Daily   --> Initialised DataFrame in loop
               raw_data={'biais_obs':[''],'lon':[''],'lat':['']}
               df_diff_day = pd.DataFrame(raw_data, columns = ['biais_obs','lon','lat'])
                
                
                 
            #cterm recuperer heure pour bouees
               cterm=date_reseau.hour
               cterm_mer=date_valid.hour
               print('*************date_reseau',date_reseau,'heure_reseau',cterm,'ech',cterm_for[h],'date_valid_mer',date_valid,'heure_valid_mer',cterm_mer)
               #Date xp
               date_2 = date_reseau.strftime("%d/%m/%y %H:%M")
               ##Date bouees
               date_2_bouees = date_valid.strftime("%d/%m/%y %H:%M")
                
                
            
               #Write date format repertory name for buyos
               #Date - Repertory - Buyos
               date_buyos_dir = '%4.4i%2.2i%2.2iT%2.2i00A'%(date_valid.year,\
               date_valid.month,\
               date_valid.day,\
               date_valid.hour)
               print('Date fic NC',date_buyos_dir)
                
               #Read NetCDF
               #chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(expe[e])+'/'
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
                 
                    # ind_flag = np.where((odb_key_val_dic['col_10'] ==1))
            
                    # Flag zone --> Zone (Do dict) N.atl
                    ind_flag = np.where((odb_key_val_dic['col_1'] > zones[key][0]) &
                            (odb_key_val_dic['col_1'] < zones[key][1]) &
                            (odb_key_val_dic['col_2'] > zones[key][2]) &
                            (odb_key_val_dic['col_2'] < zones[key][3]) &
                            (odb_key_val_dic['col_10'] ==1))

                    ind_flag_obs =  np.where((odb_key_val_dic['col_1'] > zones[key][0]) &
                            (odb_key_val_dic['col_1'] < zones[key][1]) &
                            (odb_key_val_dic['col_2'] > zones[key][2]) &
                            (odb_key_val_dic['col_2'] < zones[key][3]))
        
        
        
                    #an-depar
                    data2plotc = -odb_key_val_dic['col_9'][ind_flag]    ## 9 an-depar (bouees-analyse), 8 fg-depar
                    lon        = odb_key_val_dic['col_1'][ind_flag]
                    lat        = odb_key_val_dic['col_2'][ind_flag]
                    id_an      = odb_key_val_dic['col_5'][ind_flag]
                    
                    
                 
                    #obs value datum flag=1
                    data2plota   = odb_key_val_dic['col_7'][ind_flag]## 7 obsvalue
                    lon_obs      = odb_key_val_dic['col_1'][ind_flag]
                    lat_obs      = odb_key_val_dic['col_2'][ind_flag]
                    id_obs       = odb_key_val_dic['col_5'][ind_flag]
                    date         = odb_key_val_dic['col_3'][ind_flag]      # in secondes
                    time         = odb_key_val_dic['col_4'][ind_flag]      # in secondes
                
                    #obs value total 
                    data2plota_tot   = odb_key_val_dic['col_7'][ind_flag_obs]## 7 obsvalue
                    lon_obs_tot      = odb_key_val_dic['col_1'][ind_flag_obs]
                    lat_obs_tot      = odb_key_val_dic['col_2'][ind_flag_obs]
                    id_obs_tot       = odb_key_val_dic['col_5'][ind_flag_obs]
                    date_tot         = odb_key_val_dic['col_3'][ind_flag_obs]      # in secondes
                    time_tot         = odb_key_val_dic['col_4'][ind_flag_obs]      # in secondes


                           
                    
                    mask = (lon_obs >zones[key][0]) & (lon_obs < zones[key][1]) & (lat_obs > zones[key][2]) & (lat_obs < zones[key][3])
                    
                    #DF buyos data from nc files
                    a=pd.DataFrame(lon_obs)
                    b=pd.DataFrame(lat_obs)
                    c=pd.DataFrame(data2plota-273.15)
                    date=pd.DataFrame(time)                             # secondes
                    temps=pd.DataFrame(time)
                    
                    
                    
                    #DATAFRAME BUYOS -- meme bouées pas même heure
                    df_bouees=pd.concat([a,b,c,date,temps],axis=1)
                    df_bouees.columns =['lon', 'lat', 'obs','date','time']
                    obs_sst=data2plota
                    
                    #Buyos numbers
                    Nb_obs_assim=len(data2plota)
                    Nb_obs=len(data2plota_tot)
                    Nb_diff=Nb_obs-Nb_obs_assim  #Diff between nubers of obs tot and numbers of obs assim
                     
                    
                    #Read Xp data
                    # =============================================================================
                    r_for=get_file(expe[e],date_reseau,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for[h],geometry,formats,fill_for)
                    X=r_for.geometry.get_lonlat_grid()
                    lon_expe=X[0]
                    lat_expe=X[1]
                    x=r_for.readfield(name[0])
                    
                     
                    #Loop on dataframe of obs
                    # =============================================================================
                    for i in range(Nb_obs_assim):
                        SST=x.extract_point(df_bouees.lon[i],df_bouees.lat[i])   #Extract data on lon/lat
                        # print(SST)
                        biais=SST.data-df_bouees.obs[i]
                        
                        new_row1 = {'biais_obs':biais,'lon':df_bouees.lon[i], 'lat':df_bouees.lat[i]}
                        df_diff_day.loc[len(df_diff_day)] = new_row1
                    
                    # delete firt row
                    df_diff_day.drop(index=df_diff_day.index[0], axis=0, inplace=True)
                    
                    df_diff_day.drop(df_diff_day.loc[df_diff_day['biais_obs']>1.5e+100].index, inplace=True) # delete Nan values     
                   
                    # Compute mean diff and std over wone, eache day 
                    biais_obs_mean=df_diff_day['biais_obs'].mean()
                    biais_obs_std=df_diff_day['biais_obs'].std()
                    
                    new_row = {'date_valid':date_2_bouees,'date_xp':date_2, 'expe':expe[e], 'var':name[0],'level':level[0], 'ech':cterm_for[h],'biais_obs': biais_obs_mean,'eqm': biais_obs_std,'Nb_bouees_tot':Nb_obs ,'Nb_bouees_an':Nb_obs_assim,'diff_Nb':Nb_diff}
                    print('key2',key)
                    df_bouees_final[key].loc[len(df_bouees_final[key])] = new_row

        # r_for.close()
            
## Enregistre df_separe en fct de expe et zones 
    for zon in zones.keys():
        print('son',zon)
        print('data_zon[zon]', df_bouees_final[zon])
        df_bouees_final[zon].to_csv(pathres+str(expe[e])+'_'+str(zon)+'_prevision_bouees_flag_ech_maskcorr_fevrier.csv')

