#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 18:36:06 2023

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


# name=['SFX.SIC','SFX.SST','SFX.ICEHSI_1','SFX.ICETSF_1']
#name=['SFX.SIC']
#name=['SFX.SST']
# MINMAX=[-0.5,0.5]
# MINMAX=[0,1]
# MINMAX=[270,280] #rajout
# MINMAX=[[0,1],[270,280],[0,4.],[-40,1]]
#name = ['SFX.TEMP_OC21']  # 'SFX.TEMP_OC2','SFX.TEMP_OC3','SFX.TEMP_OC4','SFX.TEMP_OC5','SFX.TEMP_OC6','SFX.TEMP_OC7','SFX.TEMP_OC8','SFX.TEMP_OC9','SFX.TEMP_OC10','SFX.TEMP_OC11','SFX.TEMP_OC12']
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
cblock_oc_drown = 'c933_1'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
cmodel = 'surfex'
cmodel_oc = 'arpege'

# chemin='/d0/images/ormieresl/'
chemin = '/home/ormieresl/Scripts/'
chemin2 = '/home/ormieresl/Scripts/cartes_diffoc/'
chemin_sst='/home/ormieresl/'

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
# expe,date,cblock_oc,ckind_oc,cfields  --> dqns get file 2 arguments voir pour def les 5 arguments

r_cmo = get_file(expe, date, cblock, ckind, cfields, cmodel)
r_cmo_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel)


r = get_file(expe, date, cblock_oc, ckind_oc, cfields_oc, cmodel_oc)
r_oc = get_file(expe_test, date, cblock_oc_drown, ckind_oc, cfields_oc, cmodel_oc)
# for i in range(len(name)):
#     print(i, name[i])

f_cmo = r_cmo.readfield(name[0])



f = r.readfield(name_mer[0])
f_drown = r_oc.readfield(name_mer[0])





X = r.geometry.get_lonlat_grid()
lon=X[0]
lat=X[1]
 

mask=(f.data<1e22)
mask_drown=(f_drown.data<100000000000000000000000000000)
## Mask def
mask_def=(f.data > 200) & (f.data < 400)
mask_def_drown=(f_drown.data > 200) & (f_drown.data < 400)

# Maks undef 
mask_undef=(f.data<200)
mask_undef_drown=(f_drown.data<200)

print('total sans drown',len(f.data[mask]))
print('total avec drown',len(f_drown.data[mask_drown]))
# total sans drown 4657665
# total avec drown 4666624


print('def',len(f.data[mask_def]))
print('def drown',len(f_drown.data[mask_def_drown]))
## Les deux fichier avec et sans drown ocean ont même nb total de def
# def 2476801
# def drown 2485760

print('undef',len(f.data[mask_undef]))
print('undef drown',len(f_drown.data[mask_undef_drown]))
## Les deux fichier avec et sans drown ocean ont même nb total de undef
# undef 2180864
# undef drown 2180864

## Pour retenir que là où race.sea != 0
#f_diff = r_cmo.readfield(name[0])-(r_cmo_ref.readfield(name[0])) 
f_diff = r_cmo.readfield(name[0]) 


print('Def tot Sans drown',len(f_diff.data[mask_def]))
print('Def tot Avec drown',len(f_diff.data[mask_def_drown]))

print('Undef tot Sans drown',len(f_diff.data[mask_undef]))
print('Undef tot Avec drown',len(f_diff.data[mask_undef_drown]))


diff_def_nb=len(f_diff.data[mask_def]) - len(f_diff.data[mask_def_drown])
# Diff de points = 8959, expe_test: GMPE
# Diff de points = 0, expe_test: GM6E
print('diff de def',diff_def_nb)



### Points libres: Undef du à Mercator Undef
# print(min(f_diff.data.any()))

## Fichier arpege

mask_mer=(f_diff.data>-15) & (f_diff.data<50) ## 4679952 def
#mask_mer=(f_diff.data>-2)                     ## 6478817 total
#mask_mer=(f_diff.data<100000000000000000000000000000) ## 4679952

mask_mer_paslibre=(f_diff.data>-15) & (f_diff.data<50) & (f_drown.data > 200) & (f_drown.data < 400)

print('tous pts definis de f_diff',len(f_diff.data[mask_mer]))
print('tous pts definis de f_diff ou Mercator est definir',len(f_diff.data[mask_mer_paslibre]))
# print(len(f_diff.data[mask_def_drown]))
# print(len((f_diff.data[mask_def])))
# print(len(f_diff.data[mask_undef_drown]))


## Cartes

# print(np.where(f_diff.data>999.0))
# print(f_diff.data[0])
# print(np.shape(f_diff.data))
# print(type(f_diff.data))


# f_diff.data.astype(int)
# print(f_diff.data)


##  Cartes sans forcer les valeurs def et non def

print('points libres',len(f_diff.data[mask_undef]))
print('points non libres',len(f_diff.data[mask_def]))

chemin_undef='/home/ormieresl/Scripts/undef/'
print(chemin_undef)

crs = ccrs.PlateCarree()
# crs=ccrs.SouthPolarStereo(central_longitude=0.0, globe=None)
fig, ax = f_diff.cartoplot(
projection=crs, plot_method='scatter', colormap='seismic')
#ax.set_extent([60, 85, 0, 30], ccrs.PlateCarree())
ax.title.set_text(name[0]+'_'+str(expe)+'-'+ str(xpref)+
  '_Mask_sansdrown'+str(date)[0:10]+'-'+str(date)[11:13])
#ax.set_extent([-110, -50, 40, 70], ccrs.PlateCarree())
ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
fig.savefig(chemin_undef+name[0]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date)
   [0:10]+'-'+str(date)[11:13]+'-GLOB_sans_drowning-'+cblock+'-'+ckind+str(term)+'.png')





### Cartes undef sans drown

f_diff.data[mask_undef]=500
f_diff.data[mask_def]=-999

## 
print('points libres',len(f_diff.data[mask_undef]))
print('points non libres',len(f_diff.data[mask_def]))

crs = ccrs.PlateCarree()
# crs=ccrs.SouthPolarStereo(central_longitude=0.0, globe=None)
fig, ax = f_diff.cartoplot(
projection=crs, plot_method='scatter', colormap='seismic')
#ax.set_extent([60, 85, 0, 30], ccrs.PlateCarree())
ax.title.set_text(name[0]+'_'+str(expe)+'-'+ str(xpref)+
  '_Mask_sansdrown'+str(date)[0:10]+'-'+str(date)[11:13])
#ax.set_extent([-110, -50, 40, 70], ccrs.PlateCarree())
ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
fig.savefig(chemin_undef+name[0]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date)
   [0:10]+'-'+str(date)[11:13]+'-GLOB_sans_drowning-'+cblock+'-'+ckind+str(term)+'.png')


### Cartes undef avec drown

f_diff.data[mask_undef_drown]=500
f_diff.data[mask_def_drown]=-999

print('points libres avec drown',len(f_diff.data[mask_undef_drown]))
print('points non libres sans drown',len(f_diff.data[mask_def_drown]))


crs = ccrs.PlateCarree()
# crs=ccrs.SouthPolarStereo(central_longitude=0.0, globe=None)
fig, ax = f_diff.cartoplot(
projection=crs, plot_method='scatter', colormap='seismic')
#ax.set_extent([60, 85, 0, 30], ccrs.PlateCarree())
ax.title.set_text(name[0]+'_'+str(expe) +'-' + str(xpref)+
  '_Mask_avecdrown'+str(date)[0:10]+'-'+str(date)[11:13])
#ax.set_extent([-110, -50, 40, 70], ccrs.PlateCarree())
ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
fig.savefig(chemin_undef+name[0]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date)
   [0:10]+'-'+str(date)[11:13]+'-GMPE_drowning-'+cblock+'-'+ckind+str(term)+'.png')









# print(np.where(f_diff.data>2))
# print(len(f_diff.data==2))

# print(f_diff.data[0])
# print(np.shape(f_diff.data))
# print(f_diff.data)
# print(type(f_diff.data))
# # f_diff.data.getdata(-2.0)


# f_diff.data.astype(int)
# print(f_diff.data)