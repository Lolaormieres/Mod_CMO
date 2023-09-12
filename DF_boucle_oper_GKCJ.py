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
# name_mer = ['SURFSST.CLIM.']
# name = ['SFX.TEMP_OC1', 'SFX.TEMP_OC2', 'SFX.TEMP_OC3', 'SFX.TEMP_OC4', 'SFX.TEMP_OC5', 'SFX.TEMP_OC6', 'SFX.TEMP_OC7', 'SFX.TEMP_OC8', 'SFX.TEMP_OC9', 'SFX.TEMP_OC10', 'SFX.TEMP_OC11', 'SFX.TEMP_OC12']#,
#         #'SFX.TEMP_OC13','SFX.TEMP_OC14', 'SFX.TEMP_OC15', 'SFX.TEMP_OC16', 'SFX.TEMP_OC17', 'SFX.TEMP_OC18', 'SFX.TEMP_OC19', 'SFX.TEMP_OC20', 'SFX.TEMP_OC21', 'SFX.TEMP_OC22', 'SFX.TEMP_OC23', 'SFX.TEMP_OC24', 'SFX.TEMP_OC25', 'SFX.TEMP_OC26']
name_sal = ['SFX.SALT_OC1', 'SFX.SALT_OC2', 'SFX.SALT_OC3', 'SFX.SALT_OC4', 'SFX.SALT_OC5', 'SFX.SALT_OC6', 'SFX.SALT_OC7', 'SFX.SALT_OC8', 'SFX.SALT_OC9', 'SFX.SALT_OC10', 'SFX.SALT_OC11', 'SFX.SALT_OC12', 'SFX.SALT_OC13']
#             'SFX.SALT_OC14', 'SFX.SALT_OC15', 'SFX.SALT_OC16', 'SFX.SALT_OC17', 'SFX.SALT_OC18', 'SFX.SALT_OC19', 'SFX.SALT_OC20', 'SFX.SALT_OC21', 'SFX.SALT_OC22', 'SFX.SALT_OC23', 'SFX.SALT_OC24', 'SFX.SALT_OC25', 'SFX.SALT_OC26']
name_mer = ['SURF.THETAO1']#, 'SURF.THETAO2', 'SURF.THETAO3', 'SURF.THETAO4', 'SURF.THETAO5', 'SURF.THETAO6', 'SURF.THETAO7', 'SURF.THETAO8', 'SURF.THETAO9', 'SURF.THETAO10', 'SURF.THETAO11', 'SURF.THETAO12']#,
#            # 'SURF.THETAO13','SURF.THETAO14', 'SURF.THETAO15', 'SURF.THETAO16', 'SURF.THETAO17', 'SURF.THETAO18', 'SURF.THETAO19', 'SURF.THETAO20', 'SURF.THETAO21', 'SURF.THETAO22', 'SURF.THETAO23', 'SURF.THETAO24', 'SURF.THETAO25', 'SURF.THETAO26']
name_sal_mer = ['SURF.SALINO1']#, 'SURF.SALINO2', 'SURF.SALINO3', 'SURF.SALINO4', 'SURF.SALINO5', 'SURF.SALINO6', 'SURF.SALINO7', 'SURF.SALINO8', 'SURF.SALINO9', 'SURF.SALINO10', 'SURF.SALINO11', 'SURF.SALINO12', 'SURF.SALINO13',
#                 'SURF.SALINO14', 'SURF.SALINO15', 'SURF.SALINO16', 'SURF.SALINO17', 'SURF.SALINO18', 'SURF.SALINO19', 'SURF.SALINO20', 'SURF.SALINO21', 'SURF.SALINO22', 'SURF.SALINO23', 'SURF.SALINO24', 'SURF.SALINO25', 'SURF.SALINO26']
# name_u = ['SFX.UCUR_OC1', 'SFX.UCUR_OC2', 'SFX.UCUR_OC3', 'SFX.UCUR_OC4', 'SFX.UCUR_OC5', 'SFX.UCUR_OC6', 'SFX.UCUR_OC7', 'SFX.UCUR_OC8', 'SFX.UCUR_OC9', 'SFX.UCUR_OC10', 'SFX.UCUR_OC11', 'SFX.UCUR_OC12', 'SFX.UCUR_OC13',
#           'SFX.UCUR_OC14', 'SFX.UCUR_OC15', 'SFX.UCUR_OC16', 'SFX.UCUR_OC17', 'SFX.UCUR_OC18', 'SFX.UCUR_OC19', 'SFX.UCUR_OC20', 'SFX.UCUR_OC21', 'SFX.UCUR_OC22', 'SFX.UCUR_OC23', 'SFX.UCUR_OC24', 'SFX.UCUR_OC25', 'SFX.UCUR_OC26']
# name_v = ['SFX.VCUR_OC1', 'SFX.VCUR_OC2', 'SFX.VCUR_OC3', 'SFX.VCUR_OC4', 'SFX.VCUR_OC5', 'SFX.VCUR_OC6', 'SFX.VCUR_OC7', 'SFX.VCUR_OC8', 'SFX.VCUR_OC9', 'SFX.VCUR_OC10', 'SFX.VCUR_OC11', 'SFX.VCUR_OC12', 'SFX.VCUR_OC13',
#           'SFX.VCUR_OC14', 'SFX.VCUR_OC15', 'SFX.VCUR_OC16', 'SFX.VCUR_OC17', 'SFX.VCUR_OC18', 'SFX.VCUR_OC19', 'SFX.VCUR_OC20', 'SFX.VCUR_OC21', 'SFX.VCUR_OC22', 'SFX.VCUR_OC23', 'SFX.VCUR_OC24', 'SFX.VCUR_OC25', 'SFX.VCUR_OC26']

## Lecture des niveaux de la CMO
filename = 'level2.txt'
level = np.loadtxt(filename, delimiter=',', dtype=float)
print('level', level)
print(type(level))
print(np.shape(level))
level_12 = level[0:12]
print(level_12)
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

nb=106

# =============================================================================
# ## Name variables & dates
# =============================================================================
# term='24'
level=level_12
zone='nord_atlantique'
compar='an-mer_sst'
ech=cterm
# GKCJ (Double) V5
# GK7C (Double de nos exp) v4
geometry = 'global1798'

directory='DF_SST_REF'
expe=['GK7C','GKCJ','GN84']
# =============================================================================
# ## Lectur ficher fa
# =============================================================================

def get_file(xpid, date, rep, typfic, fields, mod, chain, ech):
    # MODEL
    fp = {'experiment': xpid, 'kind': typfic, 'block': rep, 'fields': fields,
          'date': date, 'term': ech, 'geometry': geometry, 'format': 'fa',
          'local': 'tmp2.fic', 'cutoff': chain, 'vapp': 'arpege',
          'vconf': '4dvarfr', 'model': mod, 'namespace': 'vortex.multi.fr', 'filling': 'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmp2.fic', openmode='r')
    os.remove('tmp2.fic')
    #return fcref.contents.data
    #fcref.container.filename
    return r


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
                'compar': ['']}
df = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar'])

## Remarque: Compare toujours à GL55
# =============================================================================
### Boucle sur exp, date + niveaux
#for e in range(len(expe)):
# Date
# -------------------------------
for e in range(len(expe)):
    date_init = datetime(2022, 9,1 , 0)
    # 02/09 -> 21/09 - 20 JOURS
    
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
            # print(e, expe[e])
            
            try:
                r = get_file(expe[0], date, cblock, ckind, cfields, cmodel, cut, cterm_oper)
                # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                r_mer = get_file(expe[2], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
                
                # f=r.readfield(name[i])
                # f_mer=r_mer.readfield(name_mer[i])
                # f_for=r_for.readfield(name[i])
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
                 
                 
                f_diff_mean = f_diff.data[mask].mean()
                print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                #print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[0], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm,'zone':zone,'compar':compar}
                # Use the loc method to add the new row to the DataFrame
                df.loc[len(df)] = new_row
            except:
                print('Pas de fichier')
    # =============================================================================
    # ## Sauvegarde dataframe 
    
    ## Supprimer premiere ligne duf Df
    #df.drop(index=df.index[0], axis=0, inplace=True) 
    #df2 = df.copy()
    #print(df)
df.to_csv(str(directory)+'/'+'double_'+str(compar)+'_'+str(zone)+'_.csv')

#%%
# =============================================================================
# ### Previ - Mercator 
# =============================================================================
ech='96'
compar='prev-mer_sst'  
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
df2 = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar'])

## FAire GKPH à la fin
#expe=['GK7C','OPER','GL55','GKOR','GM6E','GMOT'] ## Eventuellement boucle sur exp



#for e in range(len(expe)):
# Date
# -------------------------------
for e in range(len(expe)):
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
                r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_oper) ## Forecast
                r_mer = get_file(Mer, date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm) ## Mercator
                
                
                f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i])))
                
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                
                mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                #mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                #mask = (lat < 20) & (lat > -20) # Tropiques
                
                
                f_diff_mean = f_diff.data[mask].mean()
                print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                # Define the new row to be added
                
                new_row = {'date':date_2, 'expe':expe[0], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm_for,'zone':zone,'compar':compar}
                
                # Use the loc method to add the new row to the DataFrame
                df2.loc[len(df2)] = new_row
            except:
                print('Pas de fichier')    
    # =============================================================================
    # ## Sauvegarde dataframe 
    ## Supprimer premiere ligne duf Df
    #df2.drop(index=df.index[0], axis=0, inplace=True)
    #print(df2)
    # Df avec tous à la fin
    #df2.to_csv('DF'+'/'+'septembre_'+str(compar)+'_'+str(zone)+'_'+str(cterm_for)+'.csv')
df2.to_csv(str(directory)+'/'+'double_'+str(compar)+'_'+str(zone)+'_.csv')

#%%
# =============================================================================
# ## Name variables & dates
# =============================================================================
# term='24'
level=level_12
zone='eurat'
compar='an-mer_sst'
ech=cterm

# expe = 'GM6E'
## Definition color
#color_exp='#648FFF'  ## GL55
#color_exp='#785EF0'  ## GKOR
#color_exp='#FE6100'  ## GKPH
#color_exp='grey'  ## GM6E

# geometry='global798c22'
geometry = 'global1798'
# chemin='/d0/images/ormieresl/'
chemin = '/home/ormieresl/Scripts/Dataframe_mer/'


# #%%
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
                'compar': ['']}
df = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar'])


## Enlève GKPH va jusqu'au 27 sept, Carole prolonge OK
#expe=['GK7C','OPER','GL55','GKOR','GM6E','GMOT','GKPH'] ## Eventuellement boucle sur exp
# expe=['GK7C','GKCJ','GN84'] 

## Remarque: Compare toujours à GL55
# =============================================================================
### Boucle sur exp, date + niveaux


#for e in range(len(expe)):
# Date
# -------------------------------
for e in range(len(expe)):
    date_init = datetime(2022, 9,1 , 0)
    # 02/09 -> 21/09 - 20 JOURS
    
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
            # print(e, expe[e])
            
            try:
                r = get_file(expe[e], date, cblock, ckind, cfields, cmodel, cut, cterm_oper)
                # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                r_mer = get_file(expe[2], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
                
                # f=r.readfield(name[i])
                # f_mer=r_mer.readfield(name_mer[i])
                # f_for=r_for.readfield(name[i])
                f_diff = r.readfield(name[i])-(r_mer.readfield(name_mer[i]))
                
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                #mask=(lon>6)&(lon<7)&(lat>59)&(lat<60) #masques sur lon et lat
                #mask = (lon < 10) & (lon > -70) & (lat <-20) & (lat >-70)
                #mask=(lon>6)&(lon<7)&(lat>59)&(lat<60) #masques sur lon et lat
                #mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                # mask = (lat < 20) & (lat > -20) # Tropiques
                 
                 
                f_diff_mean = f_diff.data[mask].mean()
                #print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                #print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm,'zone':zone,'compar':compar}
                # Use the loc method to add the new row to the DataFrame
                df.loc[len(df)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 

## Supprimer premiere ligne duf Df
#df.drop(index=df.index[0], axis=0, inplace=True) 
#df2 = df.copy()
#print(df)

 
df.to_csv(str(directory)+'/'+'double_'+str(compar)+'_'+str(zone)+'_.csv')
    

#%%
# =============================================================================
# ### Previ - Mercator 
# =============================================================================
ech='96'
compar='prev-mer_sst'  
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
df2 = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar'])

## FAire GKPH à la fin
#expe=['GK7C','OPER','GL55','GKOR','GM6E','GMOT'] ## Eventuellement boucle sur exp
#expe=['GM6E','GMOT'] 
# expe=['GK7C','GKCJ','GL55'] 
#for e in range(len(expe)):
# Date
# -------------------------------

for e in range(len(expe)):
    date_init = datetime(2022, 9,1 , 0)
    # 02/09 -> 21/09 - 20 JOURS
    
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
                r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_oper) ## Forecast
                r_mer = get_file(expe[2], date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm) ## Mercator
                
                
                f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i])))
                
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                
                #mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                # mask = (lat < 20) & (lat > -20) # Tropiques
                
                
                f_diff_mean = f_diff.data[mask].mean()
                print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[0], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm_for,'zone':zone,'compar':compar}
                
                # Use the loc method to add the new row to the DataFrame
                df2.loc[len(df2)] = new_row
            
            except:
                print('Pas de fichier')
        
    # =============================================================================
    # ## Sauvegarde dataframe 
    ## Supprimer premiere ligne duf Df
    #df2.drop(index=df.index[0], axis=0, inplace=True)
#print(df2)
df2.to_csv(str(directory)+'/'+'double_'+str(compar)+'_'+str(zone)+'_.csv')
           
           
%%
=============================================================================
## Name variables & dates
=============================================================================
term='24'
level=level_12
zone='tropiques'
compar='an-mer_sst'
ech=cterm

# expe = 'GM6E'
## Definition color
#color_exp='#648FFF'  ## GL55
#color_exp='#785EF0'  ## GKOR
#color_exp='#FE6100'  ## GKPH
#color_exp='grey'  ## GM6E

# geometry='global798c22'
geometry = 'global1798'
# chemin='/d0/images/ormieresl/'
chemin = '/home/ormieresl/Scripts/Dataframe_mer/'


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
                'compar': ['']}
df = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar'])


## Enlève GKPH va jusqu'au 27 sept, Carole prolonge OK
#expe=['GK7C','OPER','GL55','GKOR','GM6E','GMOT','GKPH'] ## Eventuellement boucle sur exp
# expe=['GK7C','GKCJ','GL55']

## Remarque: Compare toujours à GL55
# =============================================================================
### Boucle sur exp, date + niveaux


#for e in range(len(expe)):
# Date
# -------------------------------
for e in range(len(expe)):
    date_init = datetime(2022, 9,1 , 0)
    # 02/09 -> 21/09 - 20 JOURS
    
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
            # print(e, expe[e])
            
            try:
                r = get_file(expe[0], date, cblock, ckind, cfields, cmodel, cut, cterm_oper)
                # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
                r_mer = get_file(expe[2], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
                
                # f=r.readfield(name[i])
                # f_mer=r_mer.readfield(name_mer[i])
                # f_for=r_for.readfield(name[i])
                f_diff = r.readfield(name[i])-(r_mer.readfield(name_mer[i]))
                
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                #mask=(lon>6)&(lon<7)&(lat>59)&(lat<60) #masques sur lon et lat
                #mask = (lon < 10) & (lon > -70) & (lat <-20) & (lat >-70)
                #mask=(lon>6)&(lon<7)&(lat>59)&(lat<60) #masques sur lon et lat
                #mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                # mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                mask = (lat < 20) & (lat > -20) # Tropiques
                 
                 
                f_diff_mean = f_diff.data[mask].mean()
                #print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                #print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm,'zone':zone,'compar':compar}
                # Use the loc method to add the new row to the DataFrame
                df.loc[len(df)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 

## Supprimer premiere ligne duf Df
#df.drop(index=df.index[0], axis=0, inplace=True) 
#df2 = df.copy()
#print(df)

df.to_csv(str(directory)+'/'+'double_'+str(compar)+'_'+str(zone)+'_.csv')
%%
# =============================================================================
# ### Previ - Mercator 
# =============================================================================
ech='96'
compar='prev-mer_sst'  
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
df2 = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','eqm','ech','zone','compar'])

## FAire GKPH à la fin
#expe=['GK7C','OPER','GL55','GKOR','GM6E','GMOT'] ## Eventuellement boucle sur exp
# expe=['GK7C','GKCJ','GL55'] 

# for e in range(len(expe)):
# Date
# -------------------------------

for e in range(len(expe)):
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
                r_for=get_file(expe[0],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_oper) ## Forecast
                r_mer = get_file(expe[2], date_valid, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm) ## Mercator
                
                
                f_diff = r_for.readfield(name[i])-((r_mer.readfield(name_mer[i])))
                
                X = r_mer.geometry.get_lonlat_grid()  # lon et lat
                lon = X[0]
                lat = X[1]
                long = 0
                lati = 20
                
                #mask = (lon < 0) & (lon > -80) & (lat > 20) & (lat < 70)  # N atlantique
                # mask = (lon >-35) & (lon < 45) & (lat > 20) & (lat < 72)  # Eurat
                mask = (lat < 20) & (lat > -20) # Tropiques
                
                
                f_diff_mean = f_diff.data[mask].mean()
                print('mean',f_diff_mean)
                f_diff_std = f_diff.data[mask].std()
                print('std',f_diff_std)
                # Date
                date_2 = date.strftime("%d/%m/%y")
                # Define the new row to be added
                new_row = {'date':date_2, 'expe':expe[e], 'var':name[i],'level':level[i], 'biais':f_diff_mean,'eqm':f_diff_std, 'ech':cterm_for,'zone':zone,'compar':compar}
                
                # Use the loc method to add the new row to the DataFrame
                df2.loc[len(df2)] = new_row
            except:
                print('Pas de fichier')
# =============================================================================
# ## Sauvegarde dataframe 
## Supprimer premiere ligne duf Df
#df2.drop(index=df.index[0], axis=0, inplace=True)
#print(df2)
# Df avec tous à la fin
df2.to_csv(str(directory)+'/'+'double_'+str(compar)+'_'+str(zone)+'_.csv')