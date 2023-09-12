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


# name=['SFX.SIC','SFX.SST','SFX.ICEHSI_1','SFX.ICETSF_1']
# MINMAX=[0,1]
# MINMAX=[270,280] #rajout
# MINMAX=[[0,1],[270,280],[0,4.],[-40,1]]
# MINMAX=[[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]]

name = ['SFX.SST']
#name_mer = ['SURFSST.CLIM.']
# name = ['SFX.TEMP_OC1', 'SFX.TEMP_OC2', 'SFX.TEMP_OC3', 'SFX.TEMP_OC4', 'SFX.TEMP_OC5', 'SFX.TEMP_OC6', 'SFX.TEMP_OC7', 'SFX.TEMP_OC8', 'SFX.TEMP_OC9', 'SFX.TEMP_OC10', 'SFX.TEMP_OC11', 'SFX.TEMP_OC12']#,
#         #'SFX.TEMP_OC13','SFX.TEMP_OC14', 'SFX.TEMP_OC15', 'SFX.TEMP_OC16', 'SFX.TEMP_OC17', 'SFX.TEMP_OC18', 'SFX.TEMP_OC19', 'SFX.TEMP_OC20', 'SFX.TEMP_OC21', 'SFX.TEMP_OC22', 'SFX.TEMP_OC23', 'SFX.TEMP_OC24', 'SFX.TEMP_OC25', 'SFX.TEMP_OC26']
# name_sal = ['SFX.SALT_OC1', 'SFX.SALT_OC2', 'SFX.SALT_OC3', 'SFX.SALT_OC4', 'SFX.SALT_OC5', 'SFX.SALT_OC6', 'SFX.SALT_OC7', 'SFX.SALT_OC8', 'SFX.SALT_OC9', 'SFX.SALT_OC10', 'SFX.SALT_OC11', 'SFX.SALT_OC12', 'SFX.SALT_OC13',
#             'SFX.SALT_OC14', 'SFX.SALT_OC15', 'SFX.SALT_OC16', 'SFX.SALT_OC17', 'SFX.SALT_OC18', 'SFX.SALT_OC19', 'SFX.SALT_OC20', 'SFX.SALT_OC21', 'SFX.SALT_OC22', 'SFX.SALT_OC23', 'SFX.SALT_OC24', 'SFX.SALT_OC25', 'SFX.SALT_OC26']
name_mer = ['SURF.THETAO1']#, 'SURF.THETAO2', 'SURF.THETAO3', 'SURF.THETAO4', 'SURF.THETAO5', 'SURF.THETAO6', 'SURF.THETAO7', 'SURF.THETAO8', 'SURF.THETAO9', 'SURF.THETAO10', 'SURF.THETAO11', 'SURF.THETAO12']#,
#            # 'SURF.THETAO13','SURF.THETAO14', 'SURF.THETAO15', 'SURF.THETAO16', 'SURF.THETAO17', 'SURF.THETAO18', 'SURF.THETAO19', 'SURF.THETAO20', 'SURF.THETAO21', 'SURF.THETAO22', 'SURF.THETAO23', 'SURF.THETAO24', 'SURF.THETAO25', 'SURF.THETAO26']
# name_sal_mer = ['SURF.SALINO1', 'SURF.SALINO2', 'SURF.SALINO3', 'SURF.SALINO4', 'SURF.SALINO5', 'SURF.SALINO6', 'SURF.SALINO7', 'SURF.SALINO8', 'SURF.SALINO9', 'SURF.SALINO10', 'SURF.SALINO11', 'SURF.SALINO12', 'SURF.SALINO13',
#                 'SURF.SALINO14', 'SURF.SALINO15', 'SURF.SALINO16', 'SURF.SALINO17', 'SURF.SALINO18', 'SURF.SALINO19', 'SURF.SALINO20', 'SURF.SALINO21', 'SURF.SALINO22', 'SURF.SALINO23', 'SURF.SALINO24', 'SURF.SALINO25', 'SURF.SALINO26']
# name_u = ['SFX.UCUR_OC1', 'SFX.UCUR_OC2', 'SFX.UCUR_OC3', 'SFX.UCUR_OC4', 'SFX.UCUR_OC5', 'SFX.UCUR_OC6', 'SFX.UCUR_OC7', 'SFX.UCUR_OC8', 'SFX.UCUR_OC9', 'SFX.UCUR_OC10', 'SFX.UCUR_OC11', 'SFX.UCUR_OC12', 'SFX.UCUR_OC13',
#           'SFX.UCUR_OC14', 'SFX.UCUR_OC15', 'SFX.UCUR_OC16', 'SFX.UCUR_OC17', 'SFX.UCUR_OC18', 'SFX.UCUR_OC19', 'SFX.UCUR_OC20', 'SFX.UCUR_OC21', 'SFX.UCUR_OC22', 'SFX.UCUR_OC23', 'SFX.UCUR_OC24', 'SFX.UCUR_OC25', 'SFX.UCUR_OC26']
# name_v = ['SFX.VCUR_OC1', 'SFX.VCUR_OC2', 'SFX.VCUR_OC3', 'SFX.VCUR_OC4', 'SFX.VCUR_OC5', 'SFX.VCUR_OC6', 'SFX.VCUR_OC7', 'SFX.VCUR_OC8', 'SFX.VCUR_OC9', 'SFX.VCUR_OC10', 'SFX.VCUR_OC11', 'SFX.VCUR_OC12', 'SFX.VCUR_OC13',
#           'SFX.VCUR_OC14', 'SFX.VCUR_OC15', 'SFX.VCUR_OC16', 'SFX.VCUR_OC17', 'SFX.VCUR_OC18', 'SFX.VCUR_OC19', 'SFX.VCUR_OC20', 'SFX.VCUR_OC21', 'SFX.VCUR_OC22', 'SFX.VCUR_OC23', 'SFX.VCUR_OC24', 'SFX.VCUR_OC25', 'SFX.VCUR_OC26']


## Lecture des niveaux de la CMO
filename = 'level2.txt'
level = np.loadtxt(filename, delimiter=',', dtype=float)
level_12 = level[0:12]
level_26 = level[0:26]


## ARP - assimilation
cfields = ''
cblock = 'surfan'
ckind = 'analysis'
cut = 'assim'
# cblock='forecast'
# ckind='historic'
## ARP - production
cblock_for = 'forecast'
ckind_for = 'historic'
cmodel = 'surfex'
cterm = '0'
cmodel_for = 'surfex'
cut_for = 'production'
cterm_for = '96'
ech=''
cterm_oper='0'
## Mercator 
cblock_oc = 'c933'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
#cfields_oc = 'sst'
cmodel_oc = 'arpege'
cut_oc = 'assim'

# GKCJ (Double) V5  
# GK7C (Double de nos exp) v4
# expe = 'GM6E'

geometry = 'global1798'
# chemin='/d0/images/ormieresl/'
chemin = '/home/ormieresl/Scripts/Dataframe_mer/'

# =============================================================================
# ## Lectur ficher fa
# =============================================================================

def get_file(xpid, date, rep, typfic, fields, mod, chain, ech):
    # MODEL
    fp = {'experiment': xpid, 'kind': typfic, 'block': rep, 'fields': fields,
          'date': date, 'term': ech, 'geometry': geometry, 'format': 'fa',
          'local': 'tmpsst_val.fic', 'cutoff': chain, 'vapp': 'arpege',
          'vconf': '4dvarfr', 'model': mod, 'namespace': 'vortex.multi.fr', 'filling': 'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpsst_val.fic', openmode='r')
    os.remove('tmpsst_val.fic')
    #return fcref.contents.data
    #fcref.container.filename
    return r



##directory
import os
Dir= "DF_SST_vf"
# Check if the directory exists
if not os.path.exists(Dir):
    # If it doesn't exist, create it
    os.makedirs(Dir)

## Nombre jours
nb=100

expe=['OPER','GL55','GKOR','GKPH','GN84','GNSR']
expe=['OPER','GL55','GKOR','GKPH','GN84','GNSR','GO4A','GNYG']
expe=['OPER','GN84','GO4A','GNYG']



# =============================================================================
# ## Name variables & dates
# =============================================================================
# term='24'
level=level_12
zone='nord_atlantique'
print(zone)
compar='an-mer_sst'
ech=cterm




#%%
# # =============================================================================
# # ### Annalyse - Mercator 
# # =============================================================================

### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': [''], 'SST_xp': [''], 'SST_mer': ['']}
df = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar','SST_xp', 'SST_mer'])



## Remarque: Compare 
# =============================================================================
### Boucle sur exp, date + niveaux
for e in range(len(expe)):
    print('zone:',zone)
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
            print(e, expe[e])
            print ('N.atl - an-mer')
            try:
                r = get_file(expe[e], date, cblock, ckind, cfields, cmodel, cut, cterm)
                # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                r_mer = get_file(expe[e], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
                
                SST_xp_tot=r.readfield(name[i]) - 273.15
                SST_mer_tot=r_mer.readfield(name_mer[i]) - 273.15
                
                f_diff = r.readfield(name[i])-(r_mer.readfield(name_mer[i]))
                            
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                
                #mask=(lon>6)&(lon<7)&(lat>59)&(lat<60) #masques sur lon et lat
                #mask = (lon < 10) & (lon > -70) & (lat <-20) & (lat >-70)
                #mask=(lon>6)&(lon<7)&(lat>59)&(lat<60) #masques sur lon et lat
                mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                #mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                #mask = (lat < 20) & (lat > -20) # Tropiques
                 
                SST_xp=SST_xp_tot.data[mask].mean()
                SST_mer=SST_mer_tot.data[mask].mean()
                
                f_diff_mean = f_diff.data[mask].mean()
                f_diff_std = f_diff.data[mask].std()
                
                # Date
                date_2 = date.strftime("%d/%m/%y")
                
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm,'zone':zone,'compar':compar, 'SST_xp':SST_xp,'SST_mer':SST_mer}
                # Use the loc method to add the new row to the DataFrame
                df.loc[len(df)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
df.to_csv(Dir+'/'+str(zone)+str(compar)+'.csv')

print('TERMINE 1')

#%%
# =============================================================================
# ### Previ - Mercator 
# =============================================================================
ech='96'
compar='prev-mer_sst'  
### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': [''], 'SST_xp': [''], 'SST_mer': ['']}
df2 = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar','SST_xp', 'SST_mer'])

for e in range(len(expe)):
    print('zone:',zone)
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
            print('N.atl - prev-mer')
            try:
                #r = get_file(expe, date_valid, cblock, ckind, cfields, cmodel, cut, cterm) ## Analyse
                r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for) ## Forecast
                r_mer = get_file(expe[e], date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm) ## Mercator
                
                SST_xp_tot=r.readfield(name[i]) - 273.15
                SST_mer_tot=r_mer.readfield(name_mer[i]) - 273.15
                
                f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i])))
                
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                
                mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                #mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                #mask = (lat < 20) & (lat > -20) # Tropiques
                SST_xp=SST_xp_tot.data[mask].mean()
                SST_mer=SST_mer_tot.data[mask].mean()
                
                f_diff_mean = f_diff.data[mask].mean()
                print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm_for,'zone':zone,'compar':compar, 'SST_xp':SST_xp,'SST_mer':SST_mer}
                # Use the loc method to add the new row to the DataFrame
                df2.loc[len(df2)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
df2.to_csv(Dir+'/'+str(zone)+str(compar)+'.csv')

#%%
# =============================================================================
# ## Name variables & dates
# =============================================================================
# term='24'
level=level_12
zone='eurat'
print(zone)
compar='an-mer_sst'
ech=cterm

# # =============================================================================
# # ### Annalyse - Mercator 
# # =============================================================================

### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': [''], 'SST_xp': [''], 'SST_mer': ['']}
df = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar','SST_xp', 'SST_mer'])
# =============================================================================
### Boucle sur exp, date + niveaux
for e in range(len(expe)):
    print('zone:',zone)
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
            print(e, expe[e])
            print('N.atl - an-mer')
            try:
                r = get_file(expe[e], date, cblock, ckind, cfields, cmodel, cut, cterm)
                # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                r_mer = get_file(expe[e], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
                
                SST_xp_tot=r.readfield(name[i]) - 273.15
                SST_mer_tot=r_mer.readfield(name_mer[i]) - 273.15
                
                f_diff = r.readfield(name[i])-(r_mer.readfield(name_mer[i]))
                            
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                # mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                #mask = (lat < 20) & (lat > -20) # Tropiques
                 
                SST_xp=SST_xp_tot.data[mask].mean()
                SST_mer=SST_mer_tot.data[mask].mean()
                
                f_diff_mean = f_diff.data[mask].mean()
                f_diff_std = f_diff.data[mask].std()
 
                # Date
                date_2 = date.strftime("%d/%m/%y")
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm,'zone':zone,'compar':compar, 'SST_xp':SST_xp,'SST_mer':SST_mer}
                # Use the loc method to add the new row to the DataFrame
                df.loc[len(df)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
df.to_csv(Dir+'/'+str(zone)+str(compar)+'.csv')

print('TERMINE 1')

#%%
# =============================================================================
# ### Previ - Mercator 
# =============================================================================
ech='96'
compar='prev-mer_sst'  
### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': [''], 'SST_xp': [''], 'SST_mer': ['']}
df2 = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar','SST_xp', 'SST_mer'])

for e in range(len(expe)):
    print('zone:',zone)
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
            print('N.atl - prev-mer')
            try:
                #r = get_file(expe, date_valid, cblock, ckind, cfields, cmodel, cut, cterm) ## Analyse
                r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for) ## Forecast
                r_mer = get_file(expe[e], date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm) ## Mercator
                  
                SST_xp_tot=r.readfield(name[i]) - 273.15
                SST_mer_tot=r_mer.readfield(name_mer[i]) - 273.15
                
                f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i])))
                
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                
                # mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                #mask = (lat < 20) & (lat > -20) # Tropiques
                SST_xp=SST_xp_tot.data[mask].mean()
                SST_mer=SST_mer_tot.data[mask].mean()
                
                f_diff_mean = f_diff.data[mask].mean()
                print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm_for,'zone':zone,'compar':compar, 'SST_xp':SST_xp,'SST_mer':SST_mer}
                
                # Use the loc method to add the new row to the DataFrame
                df2.loc[len(df2)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
df2.to_csv(Dir+'/'+str(zone)+str(compar)+'.csv')

#%%
# =============================================================================
# ## Name variables & dates
# =============================================================================
# term='24'
level=level_12
zone='tropiques'
print(zone)
compar='an-mer_sst'
ech=cterm

# # =============================================================================
# # ### Annalyse - Mercator 
# # =============================================================================

### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': [''], 'SST_xp': [''], 'SST_mer': ['']}
df = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar','SST_xp', 'SST_mer'])
# =============================================================================
### Boucle sur exp, date + niveaux
for e in range(len(expe)):
    print('zone:',zone)
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
            print(e, expe[e])
            print('N.atl - an-mer')
            try:
                r = get_file(expe[e], date, cblock, ckind, cfields, cmodel, cut, cterm)
                # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                r_mer = get_file(expe[e], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
                
                SST_xp_tot=r.readfield(name[i]) - 273.15
                SST_mer_tot=r_mer.readfield(name_mer[i]) - 273.15
                
                f_diff = r.readfield(name[i])-(r_mer.readfield(name_mer[i]))
                            
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                # mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                #mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                mask = (lat < 20) & (lat > -20) # Tropiques
                 
                SST_xp=SST_xp_tot.data[mask].mean()
                SST_mer=SST_mer_tot.data[mask].mean()
                
                f_diff_mean = f_diff.data[mask].mean()
                f_diff_std = f_diff.data[mask].std()
                # Date
                date_2 = date.strftime("%d/%m/%y")
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm,'zone':zone,'compar':compar, 'SST_xp':SST_xp,'SST_mer':SST_mer}
                # Use the loc method to add the new row to the DataFrame
                df.loc[len(df)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
df.to_csv(Dir+'/'+str(zone)+str(compar)+'.csv')

print('TERMINE 1')
# =============================================================================
# ### Previ - Mercator 
# =============================================================================
ech='96'
compar='prev-mer_sst'  
### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'eqm': [''],
                'ech': [''],
                'zone': [''],
                'compar': [''], 'SST_xp': [''], 'SST_mer': ['']}
df2 = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar','SST_xp', 'SST_mer'])
for e in range(len(expe)):    
    print('zone:',zone)
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
            print('N.atl - prev-mer')
            try:
                #r = get_file(expe, date_valid, cblock, ckind, cfields, cmodel, cut, cterm) ## Analyse
                r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for) ## Forecast
                r_mer = get_file(expe[e], date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm) ## Mercator
                  
                SST_xp_tot=r.readfield(name[i]) - 273.15
                SST_mer_tot=r_mer.readfield(name_mer[i])  - 273.15
                
                f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i])))
                
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                # mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                #mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                mask = (lat < 20) & (lat > -20) # Tropiques
                SST_xp=SST_xp_tot.data[mask].mean()
                SST_mer=SST_mer_tot.data[mask].mean()
                
                f_diff_mean = f_diff.data[mask].mean()
                print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm_for,'zone':zone,'compar':compar, 'SST_xp':SST_xp,'SST_mer':SST_mer}
                # Use the loc method to add the new row to the DataFrame
                df2.loc[len(df2)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
df2.to_csv(Dir+'/'+str(zone)+str(compar)+'.csv')
