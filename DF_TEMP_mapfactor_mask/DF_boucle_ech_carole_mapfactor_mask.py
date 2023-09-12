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
pathres = '/home/ormieresl/Routines/DF_TEMP_mapfactor_mask/'
# experience à evaluer
#expe=['GKOR','GKPH','GN84','GO4A','GNSR']
expe=['GN84']
# expe=['GN84','GO4A'] 
# référence
#expe_mer=['GN3C','GN3C','GO48','GO48','GO48']
expe_mer=['GO48']
compar='prev-mer'

# nombre de niveaux
nb_level=12

# date de départ
date_init = datetime(2023,2,1 , 0)
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
#----------------------------------------------------------------------------------------#
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


#----------------------------------------------------------------------------------------#
def sgmeand(dta,mpf) :
  ompf = 1./mpf   # Inverse map factor
  ompf_sum = np.sum(ompf) # sigma of 1/m (weight)
  dta_norm = np.ma.divide(dta,mpf) # Global map factor normalisation
  dta_sgmean = np.sum(dta_norm) / ompf_sum # Summation and normalisation by weight
  return dta_sgmean
#----------------------------------------------------------------------------------------#
def sgmeandm(dta,mpf,mask=None) :
  if mask is not None :
     mask_init = np.ma.getmask(dta) 
     mask_inv  = np.logical_not(mask)
     mask_loc = np.ma.mask_or(mask_inv,mask_init)  
     dta_loc = np.ma.array(dta,mask=mask_loc)
     mpf_loc = np.ma.array(mpf,mask=mask_loc)
     return sgmeand(dta_loc,mpf_loc)
  else :   
     return sgmeand(dta,mpf)
#----------------------------------------------------------------------------------------#


#%%

# =============================================================================
# ## Name variables & dates
# =============================================================================
#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80, 0, 20, 70], "eurat": [-35, 45, 20, 72], "tropiques": [-180, 180, -
    20, 20], "hn20": [-180, 180, 20, 90], "hs20": [-180, 180, -90, -20], "glob": [-180, 180, -90, 90],"med": [-3,16,30,44]}


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
                    r_mer = get_file(expe_mer[e], date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm_mer) ## Mercator
                    X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                    lon = X[0]
                    lat = X[1]

                    for i in range(len(name)):
                        print(i, name[i])
                        print(level[i])
                        f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i]))-273.15)
                        temp_exp = r_for.readfield(name[i])
                        temp_merc = r_mer.readfield(name_mer[i])-273.15
                        
                        #get map factor of stretch grid
                        map_factor = (f_diff.geometry.map_factor_field()).getdata()
                        for key in zones.keys():

                          maskdm = (lon >zones[key][0]) & (lon < zones[key][1]) & (lat > zones[key][2]) & (lat < zones[key][3])  
                          
                          #With map factor of the stretch grid
                          f_diff_mean = sgmeandm(f_diff.data,map_factor,mask=maskdm)
                          
                          #std
                          ecart         = np.square(f_diff.data-f_diff_mean)
                          variance      = sgmeandm(ecart,map_factor,mask=maskdm)
                          f_diff_std    = np.sqrt(variance)
                          
                          #temp
                          temp_exp_mean  = sgmeandm(temp_exp.data,map_factor,maskdm)
                          temp_merc_mean = sgmeandm(temp_merc.data,map_factor,maskdm)
                         
                          
                          #Without map factor of the stretch grid
                          #f_diff_mean = f_diff.data[mask].mean()
                          #f_diff_std = f_diff.data[mask].std()
                          #temp_exp_mean  = temp_exp.data[mask].mean() # relancer sans ecraser la val avec GN84 & GOJQ & GO4A
                          #temp_merc_mean = temp_merc.data[mask].mean()
                        
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
        data_zon[zon].to_csv(pathres+str(expe[e])+'_'+str(zon)+'_ech_fevrier_mars'+'_'+str(compar)+'.csv')

#%%
