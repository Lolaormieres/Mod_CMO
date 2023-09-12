#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 16:09:03 2023

@author: ormieresl
"""

#! /usr/bin/python

# coding: utf-8

# import bronx
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
pathres = "/home/labadie/SCORES_CMO/scorestot/"

# experience à evaluer
expe = ['GL55', 'GN84','GO4A','GNYG']#OPER??
expe_mer = ['GL55', 'GN84','GO4A','GNYG']

# référence

# expe_mer=['GO48']
compar = 'prev-mer'+expe_mer[0]

# directory
Dir = "DF_temp_vf"
# Check if the directory exists
if not os.path.exists(Dir):
    # If it doesn't exist, create it
    os.makedirs(Dir)

# date de depart
date_init = datetime(2022, 9, 1, 0)
# nombre de jours
nb =160 # au lieur 160
date = date_init

# échéances
ech_fin = 102
ech_ini = 0
ech_step = 6

# nombre de niveaux
nb_level = 12

# lecture des variables
name_sst = ['SFX.SST']
name = ['SFX.TEMP_OC'+str(i) for i in range(1, nb_level+1)]
name_sal = ['SFX.SALT_OC'+str(i) for i in range(1, nb_level+1)]
name_mer = ['SURF.THETAO'+str(i) for i in range(1, nb_level+1)]
name_sal_mer = ['SURF.SALINO'+str(i) for i in range(1, nb_level+1)]
name_u = ['SFX.UCUR_OC'+str(i) for i in range(1, nb_level+1)]
name_v = ['SFX.VCUR_OC'+str(i) for i in range(1, nb_level+1)]

# Lecture des niveaux de la CMO
filename = 'level2.txt'
levelfic = np.loadtxt(filename, delimiter=',', dtype=float)
print('level', levelfic)
print(type(levelfic))
print(np.shape(levelfic))
level = levelfic[0:nb_level]

########################################################
# definition des champs a extraire
########################################################
# ARP - assimilation
cfields = ''
cblock = 'surfan'
ckind = 'analysis'
cut = 'assim'
cterm = '0'
# cblock='forecast'
# ckind='historic'

# ARP - production - prevision
cblock_for = 'forecast'
ckind_for = 'historic'
cmodel = 'surfex'
cmodel_for = 'surfex'
cut_for = 'production'
cterm_for = '96'
# Mercator
cblock_oc = 'c933'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
# cfields_oc = 'sst'
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
          'local': 'tmptemp_expe.fic', 'cutoff': chain, 'vapp': 'arpege',
          'vconf': '4dvarfr', 'model': mod, 'namespace': 'vortex.multi.fr', 'filling': 'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmptemp_expe.fic', openmode='r')
    os.remove('tmptemp_expe.fic')
    # return fcref.contents.data
    # fcref.container.filename
    return r

# GKCJ (Double de nos exp) v5
# GK7C (Double de nos exp) v4
# expe=['OPER','GL55','GKOR','GKPH','GN84','GNSR','GO4A','GNYG']
# expe=['OPER','GN84','GO4A','GNYG']

# %%



# =============================================================================
# ## Name variables & dates
# =============================================================================
# domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80, 0, 20, 70], "eurat": [-35, 45, 20, 72], "tropiques": [0, 360, -
    20, 20], "hn20": [0, 360, 20, 90], "hs20": [0, 360, -90, -20], "glob": [0, 360, -90, 90]}
compar = 'an-mer_sst'
ech = cterm
dt = timedelta(days=1)


# # =============================================================================
# # ### Annalyse - Mercator
# # =============================================================================
# Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': [''], 'SST_xp': [''], 'SST_mer': ['']}

for e in range(len(expe)):
    df = {}
    for zon in zones.keys():
        df[zon] = pd.DataFrame(raw_data, columns=['date', 'expe', 'var', 'level',
                               'biais', 'eqm', 'ech', 'zone', 'compar', 'SST_xp', 'SST_mer'])
    dt = timedelta(days=1)
    date=date_init 
    period = timedelta(days=nb)
    date_end = date+period
    for t in range(nb):
            date = date + dt
            date_2 = date.strftime("%d/%m/%y ")
            print('date',date_2)
            try:
                r = get_file(expe[e], date, cblock, ckind,cfields, cmodel, cut, cterm)
                # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                r_mer = get_file(expe_mer[e], date, cblock_oc, ckind_oc, cfields_oc, cmodel_oc, cut_oc, cterm)
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
    
                lon = X[0]
                lat = X[1]
    
                for i in range(len(name)):
                    print(i, name[i])
                    # SST arp Kelvin, TempOC °C
                    SST_xp_tot=r.readfield(name[i]) 
                    SST_mer_tot=r_mer.readfield(
                        name_mer[i]) - 273.15  # Surf theta Kelvin
    
                    f_diff=r.readfield(name[i])-(r_mer.readfield(name_mer[i])-273.15)
                    for key in zones.keys():
                        #print('keyyyy',key)
    
                      mask=(lon > zones[key][0]) & (lon < zones[key][1]) & (
                            lat > zones[key][2]) & (lat < zones[key][3])
                      f_diff_std=f_diff.data[mask].std()
                      f_diff_mean=f_diff.data[mask].mean()
                           
                      SST_xp=SST_xp_tot.data[mask].mean()
                      SST_mer=SST_mer_tot.data[mask].mean()
    
                      # Define the new row to be added
                      new_row={'date': date_2, 'expe': expe[e], 'var': name[i], 'level': level[i], 'biais': f_diff_mean,
                'eqm': f_diff_std,'ech':cterm,  'zone': key, 'compar': compar, 'SST_xp': SST_xp, 'SST_mer': SST_mer}
            
                     # Use the loc method to add the new row to the DataFrame
                      df[key].loc[len(df[key])]=new_row
                r.close()
                r_mer.close()
            except:
                print('Pas de fichier')

# =============================================================================
## Enregistre df separé
for zon in zones.keys():
    print('zon',zon)
    print('df[zon]',df[zon])
    df[zon].to_csv(Dir+"/"+str(expe[e])+'_'+str(zon)+'_'+str(compar)+'.csv')

#%%

# =============================================================================
# ## Name variables & dates
# =============================================================================
# domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80, 0, 20, 70], "eurat": [-35, 45, 20, 72], "tropiques": [0, 360, -
    20, 20], "hn20": [0, 360, 20, 90], "hs20": [0, 360, -90, -20], "glob": [0, 360, -90, 90]}
compar = 'prev-mer_sst'
ech = cterm
dt = timedelta(days=1)


# # =============================================================================
# # ### Annalyse - Mercator
# # =============================================================================
# Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': [''], 'SST_xp': [''], 'SST_mer': ['']}


for e in range(len(expe)):
    df2 = {}
    for zon in zones.keys():
        df2[zon] = pd.DataFrame(raw_data, columns=['date', 'expe', 'var', 'level',
                               'biais', 'eqm', 'ech', 'zone', 'compar', 'SST_xp', 'SST_mer'])
    dt = timedelta(days=1)
    dt2 = timedelta(days=4)
    date=date_init 
    date_valid = date_init + dt2
    period = timedelta(days=nb)
    date_end = date+period
    for t in range(nb):
            date = date + dt
            date_2 = date.strftime("%d/%m/%y")
            date_valid = date + dt2
            date_valid_2 = date_valid.strftime("%d/%m/%y")
            print('date',date_2)
            print('date_valid',date_valid_2)
            try:
                #r = get_file(expe[e], date, cblock, ckind,cfields, cmodel, cut, cterm)
                r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                r_mer = get_file(expe_mer[e], date, cblock_oc, ckind_oc, cfields_oc, cmodel_oc, cut_oc, cterm)
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
    
                lon = X[0]
                lat = X[1]
    
                for i in range(len(name)):
                    print(i, name[i])
                    # SST arp Kelvin, TempOC °C
                    SST_xp_tot=r_for.readfield(name[i]) 
                    SST_mer_tot=r_mer.readfield(
                        name_mer[i]) - 273.15  # Surf theta Kelvin
    
                    f_diff=r_for.readfield(
                        name[i])-(r_mer.readfield(name_mer[i])-273.15)
                    for key in zones.keys():
                        #print('keyyyy',key)
    
                      mask=(lon > zones[key][0]) & (lon < zones[key][1]) & (
                            lat > zones[key][2]) & (lat < zones[key][3])
                      f_diff_std=f_diff.data[mask].std()
                      f_diff_mean=f_diff.data[mask].mean()
                           
                      SST_xp=SST_xp_tot.data[mask].mean()
                      SST_mer=SST_mer_tot.data[mask].mean()

    
                      # Define the new row to be added
                      new_row={'date': date_2, 'expe': expe[e], 'var': name[i], 'level': level[i], 'biais': f_diff_mean,
                'eqm': f_diff_std, 'zone': key, 'compar': compar, 'SST_xp': SST_xp, 'SST_mer': SST_mer}
            
                     # Use the loc method to add the new row to the DataFrame
                      df2[key].loc[len(df2[key])]=new_row
                r.close()
                r_mer.close()
            except:
                print('Pas de fichier')
# =============================================================================
## Enregistre df separé
    for zon in zones.keys():
        print('zon',zon)
        print('df[zon]',df2[zon])
        df2[zon].to_csv(Dir+"/"+str(expe[e])+'_'+str(zon)+'_'+str(compar)+'.csv')

