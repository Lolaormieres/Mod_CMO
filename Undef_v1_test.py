#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:00:04 2023

@author: ormieresl
"""


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
epygram.init_env()

MINMAX = [[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5],
          [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]]
# MINMAX=[[0,30],[0,30],[0,30],[0,30],[0,30],[0,30],[0,30],[0,30],[0,30],[0,30],[0,30],[0,30]]
name=['SFX.TEMP_OC1']#,'SFX.TEMP_OC2','SFX.TEMP_OC3','SFX.TEMP_OC4','SFX.TEMP_OC5','SFX.TEMP_OC6','SFX.TEMP_OC7','SFX.TEMP_OC8','SFX.TEMP_OC9','SFX.TEMP_OC10','SFX.TEMP_OC11','SFX.TEMP_OC12','SFX.TEMP_OC13','SFX.TEMP_OC14','SFX.TEMP_OC15','SFX.TEMP_OC16','SFX.TEMP_OC17','SFX.TEMP_OC18','SFX.TEMP_OC19','SFX.TEMP_OC20','SFX.TEMP_OC21','SFX.TEMP_OC22','SFX.TEMP_OC23','SFX.TEMP_OC24','SFX.TEMP_OC25','SFX.TEMP_OC26']
# ,'SURF.THETAO2','SURF.THETAO3','SURF.THETAO4','SURF.THETAO5','SURF.THETAO6','SURF.THETAO7','SURF.THETAO8','SURF.THETAO9','SURF.THETAO10','SURF.THETAO11','SURF.THETAO12']
name_mer = ['SURF.THETAO1']
name_sst=['SURFSST.CLIM.']
# name_mer=['SURF.THETAO1','SURF.THETAO2','SURF.THETAO3','SURF.THETAO4','SURF.THETAO5','SURF.THETAO6','SURF.THETAO7','SURF.THETAO8','SURF.THETAO9','SURF.THETAO10','SURF.THETAO11','SURF.THETAO12','SURF.THETAO13','SURF.THETAO14','SURF.THETAO15','SURF.THETAO16','SURF.THETAO17','SURF.THETAO18','SURF.THETAO19','SURF.THETAO20','SURF.THETAO21','SURF.THETAO22','SURF.THETAO23','SURF.THETAO24','SURF.THETAO25','SURF.THETAO26']
#MINMAX= [[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]]
# MINMAX=[[300,310]] # ,[-10,10],[-10,10]]
# MINMAX=[[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300]]

# Lecture des niveaux de la CMO
filename = 'level2.txt'
level = np.loadtxt(filename, delimiter=',', dtype=float)
# print('level',level)
level_12 = []
level_12.append(level[0:12])
level_26 = []
level_26.append(level[0:26])
level_12 = np.array(level_12)
level_12 = level_12.T
level_26 = np.array(level_26)
level_26 = level_26.T
# print('dim',level_26.ndim)
# print('level_12',level_12)
# print('level_26',level_26)


expe = 'GL55'
expe_test='GMPE'
expe2 = 'GL55'
expe3 = 'GKPH'
xpref = 'GKPH'
# geometry='global798c22'
geometry = 'global1798'
term = '0'

# arp
cfields = ''
cblock = 'surfan'
ckind = 'analysis'
# cblock='forecast'
# ckind='historic'
# ocean
cblock_oc = 'c933'
#cblock_oc_drown = 'c933_1'
cblock_oc_drown='c933_iter1_10_testinverse'
# cblock_oc_drown = 'c933'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
cmodel = 'surfex'
cmodel_oc = 'arpege'

# chemin='/d0/images/ormieresl/'
chemin = '/home/ormieresl/Scripts/'
chemin2 = '/home/ormieresl/Scripts/cartes_diffoc/'
chemin_sst='/home/ormieresl/'

#Sauvegarde images
chemin_undef = '/cnrm/recyf/Data/users/ormieresl/undef/'


## ARP
# cfields=''
# cblock='surfan'
# ckind='analysis'
# cut='assim'
# #cblock='forecast'
# #ckind='historic'
# cblock_for='forecast'
# ckind_for='historic'
# cmodel='surfex'
# cterm='0'
# cmodel_for='surfex'
# cut_for='production'
# cterm_for='48'
# cterm_p6='6'
# ## Mer
# cblock_oc='c933'
# ckind_oc='geofields'
# cfields_oc='ocean'
# #cfields_oc='sst'
# cmodel_oc='arpege'
# cut_oc='assim'
# Date

# -------------------------------
#GL55 et GKOR
#nb = 66
#date_init = datetime(2022, 7, 26, 0)
nb = 2
date_init = datetime(2022, 9, 8, 0)
date = date_init
dt = timedelta(days=1)
date_valid = date+dt
period = timedelta(days=nb)
date_end = date+period


def get_file(xpid, date, rep, typfic, fields, mod):
    # MODEL
    fp = {'experiment': xpid, 'kind': typfic, 'block': rep, 'fields': fields,
          'date': date, 'term': term, 'geometry': geometry, 'format': 'fa',
          'local': 'tmp2.fic', 'cutoff': 'assim', 'vapp': 'arpege',
          'vconf': '4dvarfr', 'model': mod, 'namespace': 'vortex.multi.fr', 'filling': 'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmp2.fic', openmode='r')
    return r



# =============================================================================
# ## Recuperer les fichiers
# =============================================================================
##fichiers cmo
r_cmo = get_file(expe, date, cblock, ckind, cfields, cmodel)

## fichier oc
r = get_file(expe, date, cblock_oc, ckind_oc, cfields_oc, cmodel_oc) # GL55
r_oc = get_file(expe_test, date, cblock_oc_drown, ckind_oc, cfields_oc, cmodel_oc) # drown GMOT


## Lon et lat
X = r_cmo.geometry.get_lonlat_grid()
lon=X[0]
lat=X[1]


mask = (lat<90) & (lat>-60) & (lon>0) & (lon<360)
## Récupere fields
f_cmo = r_cmo.readfield(name[0])
print(len(f_cmo.data[mask]))

f = r.readfield(name_mer[0])
print(len(f.data[mask]))
f_drown = r_oc.readfield(name_mer[0])
print(len(f_drown.data[mask]))

f_diff=r_cmo.readfield(name[0])-(r.readfield(name_mer[0])-273.15)

f_diff_drown=r_cmo.readfield(name[0])-(r_oc.readfield(name_mer[0])-273.15)



### test

print(np.size(np.where(f_diff_drown.data.data[mask]<200)))
print(np.size(np.where(f_diff.data.data[mask]<200)))

print(np.size(np.where(f_cmo.data.data[mask]==0.00)))
print(np.size(np.where(f_diff.data.data[mask]==0.00))) # 0 
print(np.size(np.where(f_diff_drown.data.data[mask]==0.000))) #0


print(np.size(np.where(f_diff_drown.data.data[mask]>200)))
print(np.size(np.where(f_diff.data.data[mask]>200)))


print(len(f_diff.data[mask]))
print(len(f_diff_drown.data[mask]))


print(len(f_cmo.data[mask]))
print(len(f_drown.data[mask]))

#print(np.size(np.where((f_diff_drown.data.data[mask]>200) & (f_diff_drown.data.data[mask]<400 ))))



# =============================================================================
# Lecture points libres a partir ficher masked epygram
# =============================================================================
#Lecture avec Mask H2D field, TRUE/FALSE
# Nb pt cmo def
Nb_cmo_def=len(f_cmo.data[np.invert(f_cmo.data.mask)])
print('Nb_cmo_def',Nb_cmo_def)


# Nb pt c933 def
Nb_c933_def=len(f.data[np.invert(f.data.mask)])
print('Nb_c933_def',Nb_c933_def)


# Nb pt c933 def avec drown
Nb_c933_drown_def=len(f_drown.data[np.invert(f_drown.data.mask)])
print('Nb_c933_drown_def',Nb_c933_drown_def)

## Sans drowning
Nb_pts_libres=len(f_cmo.data[np.invert(f_cmo.data.mask)])-len(f.data[np.invert(f.data.mask)])
print('Nb_pts_libres',Nb_pts_libres)


## Avec drowning
Nb_pts_libres_drown=len(f_cmo.data[np.invert(f_cmo.data.mask)])-len(f_drown.data[np.invert(f_drown.data.mask)])
print('Nb_pts_libres_drown',Nb_pts_libres_drown)



## Pourcentages pt libres
Quotas=Nb_pts_libres/Nb_cmo_def*100
print('pourcentages de pt libres avant drown:',Quotas)
## Pourcentages pt libres avec drown
Quotas=Nb_pts_libres_drown/Nb_cmo_def*100
print('pourcentages de pt libres avec drown:',Quotas)