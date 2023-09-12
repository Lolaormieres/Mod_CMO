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
expe=['GL55','GN84','GO4A','GNYG']#OPER??
expe_mer=['GL55','GN84','GO4A','GNYG']
# référence
expe_mer=['GO48']
compar='prev-mer'+expe_mer[0]


##directory
import os
Dir= "DF_vf"
# Check if the directory exists
if not os.path.exists(Dir):
    # If it doesn't exist, create it
    os.makedirs(Dir)


# date de depart
date_init = datetime(2022, 9,2 , 0)
# nombre de jours
nb=160

# échéances
ech_fin=102
ech_ini=0
ech_step=6

# nombre de niveaux
nb_level=12


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
# geometry='global798c22'
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
    #fcref.cont
#%%
# =============================================================================
# ## Name variables & dates
# =============================================================================
#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [0,360,-20,20],"hn20":[0,360,20,90],"hs20":[0,360,-90,-20],"glob":[0,360,-90,90]}
compar=an-mer_sst
level=level_12
ech=cterm

# =============================================================================
# ### Annalyse - Mercator 
# =============================================================================
### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': ['']}
for e in range(len(expe)):
    df={}
    for zon in zones.keys():
    df[zon] = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar'])
# =============================================================================
### Boucle sur exp, date + niveaux
        # Date
        # -------------------------------
        date_init = datetime(2022, 9,1 , 0)
        date = date_init
        dt = timedelta(days=1)
        date_valid = date+dt
        period = timedelta(days=nb)
        date_end = date+period    
    
        for t in range(nb):
            date = date + dt
            for i in range(len(name)):
                print(i, name[i])
                print(level[i])
                print(e, expe[e])
                try:
                    r = get_file(expe[e], date, cblock, ckind, cfields, cmodel, cut, cterm)
                    # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                    r_mer = get_file(expe[e], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
                    f_diff = r.readfield(name[i])-(r_mer.readfield(name_mer[i])-273.15)
                
                    X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                    lon = X[0]
                    lat = X[1]
                    long = 0
                    lati = 20
                
                    for key in zones.keys():
                        mask = (lon >zones[key][0]) & (lon < zones[key][1]) & (lat > zones[key][2]) & (lat < zones[key][3])
                        #mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                        #mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                        #mask = (lat < 20) & (lat > -20) # Tropiques
                        f_diff_mean = f_diff.data[mask].mean()
                        f_diff_std = f_diff.data[mask].std()
                        date_2 = date.strftime("%d/%m/%y")
                
                        # Define the new row to be added
                        new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm,'zone':key,'compar':compar}
                        # Use the loc method to add the new row to the DataFrame
                        df.loc[len(df)] = new_row
                    r.close()
                    r_mer.close()
                except:
                    print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
## Supprimer premiere ligne duf Df
    #df.to_csv('DF_octobre'+'/'+'septembre_octobre_'+str(compar)+'_'+str(zone)+'.csv')
    #df.to_csv('DF_4'+'/'+'septembre_octobre_'+'_'+str(compar)+'_'+str(zone)+'.csv')

## Enregistre df separé
    for zon in zones.keys():
        print('zon',zon)
        print('df[zon]',df[zon])
        df[zon].to_csv(Dir+"/"+str(expe[e])+'_'+str(zon)+'_'+str(compar)+'.csv')
print('TERMINE 1')     

# =============================================================================
# ### Previ - Mercator 
# =============================================================================
ech='96'
compar='prev-mer'  
zone='nord-atlantique'
## Prev - mer 
### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': ['']}

for e in range(len(expe)):
    df2={}
    for zon in zones.keys():
        df2[zon] = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar'])
        # Date
        # -------------------------------
        date_init = datetime(2022, 9,1 , 0)
        date = date_init
        dt = timedelta(days=1)
        dt2 = timedelta(days=4)
        date_valid = date+dt
        period = timedelta(days=nb)
        date_end = date+period   
        for t in range(nb):
            date = date+dt
            date_valid = date+dt2
            for i in range(len(name)):
                print(i, name[i])
                print(level[i])
                try:
                    #r = get_file(expe, date_valid, cblock, ckind, cfields, cmodel, cut, cterm) ## Analyse
                    r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for) ## Forecast
                    r_mer = get_file(expe[e], date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm) ## Mercator
                    f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i]))-273.15)
                    X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                    lon = X[0]
                    lat = X[1]
                    long = 0
                    lati = 20
                    for key in zones.keys():
                        mask = (lon >zones[key][0]) & (lon < zones[key][1]) & (lat > zones[key][2]) & (lat < zones[key][3])
 
                    #mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                    #mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                    #mask = (lat < 20) & (lat > -20) # Tropiques
                    
                        f_diff_mean = f_diff.data[mask].mean()
                        #print('mean',f_diff_mean)
                        f_diff_std = f_diff.data[mask].std()
                        #print('std',f_diff_std)
                        # Date
                        date_2 = date.strftime("%d/%m/%y")
                        # Define the new row to be added
                        new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm_for,'zone':key,'compar':compar}
                
                         # Use the loc method to add the new row to the DataFrame
                         df2[key].loc[len(df2[key])] = new_row
                    r_for.close()
                    r_mer.close()
               except:
                    print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
    #df2.to_csv('DF_1'+'/'+'septembre_octobre_GKPH_'+str(compar)+'_'+str(zone)+'.csv')
## Enregistre df separé
    for zon in zones.keys():
        print('zon',zon)
        print('df2[zon]',df2[zon])
        df2[zon].to_csv(Dir+"/"+str(expe[e])+'_'+str(zon)+'_'+str(compar)+'.csv')
                                                                                   311,80        Bot


