#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 16:09:03 2023

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
import cartopy.crs as ccrs
import cartopy
import epygram
from datetime import datetime, timedelta
import vortex
import pandas as pd
import os
import numpy.ma as ma
epygram.init_env()

# chemin de stockage des résultats
pathres="/home/labadie/SCORES_CMO/scorestot/"

# experience à evaluer
expe=['GNSR'] 
# référence
expe_mer=['GO48']
compar='prev-mer'+expe_mer[0]

# nombre de niveaux
nb_level=12

# date de départ
date_init = datetime(2022, 9,2 , 0)
# nombre de jours
nb=60

# échéances
ech_fin=102
ech_ini=0
ech_step=6

# lecture des variables
name=['SFX.TEMP_OC'+str(i) for i in range(1,nb_level+1)]
name_sal=['SFX.SALT_OC'+str(i) for i in range(1,nb_level+1)]
name_mer=['SURF.THETAO'+str(i) for i in range(1,nb_level+1)]
name_sal_mer=['SURF.SALINO'+str(i) for i in range(1,nb_level+1)]
name_u=['SFX.UCUR_OC'+str(i) for i in range(1,nb_level+1)]
name_v=['SFX.VCUR_OC'+str(i) for i in range(1,nb_level+1)]

## Lecture des niveaux de la CMO
filename = 'level2.txt'
levelfic = np.loadtxt(filename, delimiter=',', dtype=float)
print('level', levelfic)
print(type(levelfic))
print(np.shape(levelfic))
level = levelfic[0:nb_level]

########################################################
#definition des champs a extraire
########################################################
## ARP - assimilation
cfields = ''
cblock = 'surfan'
ckind = 'analysis'
cut = 'assim'
# cblock='forecast'
# ckind='historic'

## ARP - production - prevision
cblock_for = 'forecast'
ckind_for = 'historic'
cmodel = 'surfex'
cmodel_for = 'surfex'
cut_for = 'production'

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

def get_file(xpid, date, rep, typfic, fields, mod, chain, ech):
    # MODEL
    fp = {'experiment': xpid, 'kind': typfic, 'block': rep, 'fields': fields,
          'date': date, 'term': ech, 'geometry': geometry, 'format': 'fa',
          'local': 'tmpech.fic', 'cutoff': chain, 'vapp': 'arpege',
          'vconf': '4dvarfr', 'model': mod, 'namespace': 'vortex.multi.fr', 'filling': 'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpech.fic', openmode='r')
    os.remove('tmpech.fic')
    #return fcref.contents.data
    #fcref.container.filename
    return r

#%%

# =============================================================================
# ## Name variables & dates
# =============================================================================
#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [0,360,-20,20],"hn20":[0,360,20,90],"hs20":[0,360,-90,-20],"glob":[0,360,-90,90]}


# =====================date_valid.hour=========================================
# ### Previ - Mercator 
# =============================================================================
## Prev - mer 
### Ceate Dataframe
raw_data = {'date_xp': [''],
            'ech': [''],
            'date_mer': [''],
            'cterm_mer':[''],
            'expe': [''],
            'var': [''],
            'level': [''],
            'biais': [''],
            'eqm': [''],
            'zone': [''],
            'compar': [''],
            'temp_exp_mean':[''],
            'temp_merc_mean':['']}


####################
# start
####################



for e in range(len(expe)):
    data_zon={}
    for zon in zones.keys():
        data_zon[zon]=pd.DataFrame(raw_data, columns = ['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar','temp_exp_mean','temp_merc_mean'])
    # Date
    # -------------------------------
    
    dt = timedelta(days=1)
    date_reseau = date_init - dt
    
    date_valid=date_init-timedelta(hours=ech_step)
       
    period = timedelta(days=nb)
    date_end = date_reseau+period
    

    for t in range(nb):
            date_reseau = date_reseau+dt
            cterm_for= range(ech_ini,ech_fin,ech_step) 
            for h in range(len(cterm_for)):
                dt2=timedelta(hours=cterm_for[h])
                date_valid = date_reseau+dt2
                
                
                ### Cterm recuperer heure pour Mercator
                cterm=date_reseau.hour
                cterm_mer=date_valid.hour
                print('*************date_reseau',date_reseau,'heure_reseau',cterm,'ech',cterm_for[h],'date_valid_mer',date_valid,'heure_valid_mer',cterm_mer)
                #DATE XP
                date_2 = date_reseau.strftime("%d/%m/%y %H:%M")
                ## DATE MERCATOR
                date_2_mer = date_valid.strftime("%d/%m/%y %H:%M")
                
                try:
                    r_for=get_file(expe[e],date_reseau,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for[h]) ## Forecast
                    r_mer = get_file(expe_mer, date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm_mer) ## Mercator
                    X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                    lon = X[0]
                    lat = X[1]

                    for i in range(len(name)):
                        print(i, name[i])
                        print(level[i])
                        f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i]))-273.15)
                        temp_exp = r_for.readfield(name[i])
                        temp_merc = r_mer.readfield(name_mer[i])-273.15

                        for key in zones.keys():

                          mask = (lon >zones[key][0]) & (lon < zones[key][1]) & (lat > zones[key][2]) & (lat < zones[key][3])  
                    
                          f_diff_mean = f_diff.data[mask].mean()
                          f_diff_std = f_diff.data[mask].std()
                          temp_exp_mean = temp_exp.data[mask].mean()
                          temp_merc_mean = temp_merc.data[mask].mean()
                        
                         # Define the new row to be added
                          new_row = {'date_xp':date_2, 'ech':cterm_for[h],'date_mer':date_2_mer,'cterm_mer':cterm_mer, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std,'zone':key,'compar':compar,'temp_exp_mean':temp_exp_mean,'temp_merc_mean':temp_merc_mean}
                    
                         # Use the loc method to add the new row to the DataFrame
                          data_zon[key].loc[len(data_zon[key])] = new_row
                    r_for.close()
                    r_mer.close()
                    
                except:
                      print('Pas de fichier')
# =============================================================================

## Enregistre df separé
    for zon in zones.keys():
        print('zon',zon)
        print('data_zon[zon]',data_zon[zon])
        data_zon[zon].to_csv(pathres+"/"+str(expe[e])+'_'+str(zon)+'_ech_septembre_octobre'+'_'+str(compar)+'.csv')

#%%
