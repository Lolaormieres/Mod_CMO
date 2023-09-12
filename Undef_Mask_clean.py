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


expe = 'GL55'
expe_test='GMPE' #GMPE
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
cblock_oc_drown = 'c933_1' # 
# cblock_oc_drown='c933_iter1_10_testinverse' # CLOSEST - maybe same with 5 LOOP
cblock_oc_drown='c933_iter1_10' # test dont work 
cblock_oc_drown='c933_iter1_10_testITM' # NO
# cblock_oc_drown = 'c933_5' # test  dont work
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
r_cmo_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel)

## fichier oc
r = get_file(expe, date, cblock_oc, ckind_oc, cfields_oc, cmodel_oc)
r_oc = get_file(expe_test, date, cblock_oc_drown, ckind_oc, cfields_oc, cmodel_oc)


## Lon et lat
X = r_cmo.geometry.get_lonlat_grid()
lon=X[0]
lat=X[1]



## Récupere fields
f_cmo = r_cmo.readfield(name[0])
f = r.readfield(name_mer[0])
f_drown = r_oc.readfield(name_mer[0])


# =============================================================================
# # Mask, bug, du au undef=0.00
# =============================================================================

mask=(lon>0) & (lon <360) & (lat>-90) & (lat<90)

mask_antarctic=(lon>-180) & (lon <180) & (lat>-65) & (lat<90) 

# =============================================================================
# Lecture points libres a partir ficher masked epygram
# =============================================================================
#Lecture avec Mask H2D field, TRUE/FALSE
# Nb pt cmo def
## pt libres cotes vs antarciques
print('**************** pt libres total')
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


## pt libres cotes vs antarciques
print('**************** pt libres except antactic')
Nb_cmo_def=len(f_cmo.data.data[mask_antarctic][np.invert(f_cmo.data[mask_antarctic].mask)])
print('Nb_cmo_def',Nb_cmo_def)


# Nb pt c933 def
Nb_c933_def=len(f.data.data[mask_antarctic][np.invert(f.data[mask_antarctic].mask)])
print('Nb_c933_def',Nb_c933_def)


# Nb pt c933 def avec drown
Nb_c933_drown_def=len(f_drown.data.data[mask_antarctic][np.invert(f_drown.data[mask_antarctic].mask)])
print('Nb_c933_drown_def',Nb_c933_drown_def)

## Sans drowning
Nb_pts_libres=len(f_cmo.data.data[mask_antarctic][np.invert(f_cmo.data[mask_antarctic].mask)])-len(f.data.data[mask_antarctic][np.invert(f.data[mask_antarctic].mask)])
print('Nb_pts_libres',Nb_pts_libres)


## Avec drowning
Nb_pts_libres_drown=len(f_cmo.data.data[mask_antarctic][np.invert(f_cmo.data[mask_antarctic].mask)])-len(f_drown.data.data[mask_antarctic][np.invert(f_drown.data[mask_antarctic].mask)])
print('Nb_pts_libres_drown',Nb_pts_libres_drown)


## Pourcentages pt libres
Quotas=Nb_pts_libres/Nb_cmo_def*100
print('pourcentages de pt libres avant drown:',Quotas)
## Pourcentages pt libres avec drown
Quotas=Nb_pts_libres_drown/Nb_cmo_def*100
print('pourcentages de pt libres avec drown:',Quotas)




#2) lon lat


#mask_def_dom= (diff>-290) & (diff<290) & (diff!=0.00000000e+000) & (lon>-180) & (lon<180) & (lat<20) & (lat>-20)
mask_trop= (lon>-180) & (lon<180) & (lat<20) & (lat>-20)
mask_glob= (lon>-180) & (lon<180) & (lat<90) & (lat>-90)
mask_eurat= (lon>-35) & (lon<45) & (lat>20) & (lat<72)
mask_hs= (lon>-180) & (lon<180) & (lat<-20) & (lat>-90)
mask_hn= (lon>-180) & (lon<180) & (lat>20) & (lat<90)
mask_med= (lon>-3) & (lon<16) & (lat>30) & (lat<44)
mask_natl= (lon>-80) & (lon<0) & (lat>20) & (lat<70)

# TEST MERCATOR
print('****************************************************NB pts CMO')
# TEST CMO
#1. med
## DEF sea sst by domain
print('med pt:',len(f_cmo.data.data[mask_med][np.invert(f_cmo.data[mask_med].mask)]))

#Check if its still undef 0.00000000e+000 (FRAC.SEA)
a=f_cmo.data.data[mask_med][np.invert(f_cmo.data[mask_med].mask)]
print(np.where(a==0.00000000e+000))


#Check if its still undef  1.79769313e+308
b=f_cmo.data.data[mask_med][np.invert(f_cmo.data[mask_med].mask)]
print(np.where(a==1.79769313e+308))

#2. Eurat
## DEF sea sst by domain
print('eurat pt:',len(f_cmo.data.data[mask_eurat][np.invert(f_cmo.data[mask_eurat].mask)]))

#3. Natl
## DEF sea sst by domain
print('natl pt:',len(f_cmo.data.data[mask_natl][np.invert(f_cmo.data[mask_natl].mask)]))

#4. Hn
## DEF sea sst by domain
print('hn pt:',len(f_cmo.data.data[mask_hn][np.invert(f_cmo.data[mask_hn].mask)]))


#5. Hs
## DEF sea sst by domain
print('hs pt:',len(f_cmo.data.data[mask_hs][np.invert(f_cmo.data[mask_hs].mask)]))

#6. trop
## DEF sea sst by domain
print('trop pt:',len(f_cmo.data.data[mask_trop][np.invert(f_cmo.data[mask_trop].mask)]))

#7. glob
## DEF sea sst by domain
print('glob pt:',len(f_cmo.data.data[mask_glob][np.invert(f_cmo.data[mask_glob].mask)]))

# TEST MERCATOR
print('****************************************************NB pts Mercator sans drown')
#1. med
## DEF sea sst by domain
print('med pt mer:',len(f.data.data[mask_med][np.invert(f.data[mask_med].mask)]))

#Check if its still undef 0.00000000e+000 (FRAC.SEA)
a=f.data.data[mask_med][np.invert(f.data[mask_med].mask)]
print(np.where(a==0.00000000e+000))


#Check if its still undef  1.79769313e+308
b=f.data.data[mask_med][np.invert(f.data[mask_med].mask)]
print(np.where(a==1.79769313e+308))

#2. Eurat
## DEF sea sst by domain
print('eurat pt mer:',len(f.data.data[mask_eurat][np.invert(f.data[mask_eurat].mask)]))

#3. Natl
## DEF sea sst by domain
print('natl pt mer:',len(f.data.data[mask_natl][np.invert(f.data[mask_natl].mask)]))

#4. Hn
## DEF sea sst by domain
print('hn pt mer:',len(f.data.data[mask_hn][np.invert(f.data[mask_hn].mask)]))


#5. Hs
## DEF sea sst by domain
print('hs pt mer:',len(f.data.data[mask_hs][np.invert(f.data[mask_hs].mask)]))

#6. trop
## DEF sea sst by domain
print('trop pt mer:',len(f.data.data[mask_trop][np.invert(f.data[mask_trop].mask)]))

#7. glob
## DEF sea sst by domain
print('glob pt mer:',len(f.data.data[mask_glob][np.invert(f.data[mask_glob].mask)]))


# TEST MERCATOR
print('****************************************************NB pts Mercator sans drown')
#1. med
## DEF sea sst by domain
print('med pt drown:',len(f_drown.data.data[mask_med][np.invert(f_drown.data[mask_med].mask)]))

#Check if its still undef 0.00000000e+000 (FRAC.SEA)
a=f_drown.data.data[mask_med][np.invert(f_drown.data[mask_med].mask)]
print(np.where(a==0.00000000e+000))


#Check if its still undef  1.79769313e+308
b=f_drown.data.data[mask_med][np.invert(f_drown.data[mask_med].mask)]
print(np.where(a==1.79769313e+308))

#2. Eurat
## DEF sea sst by domain
print('eurat pt drown:',len(f_drown.data.data[mask_eurat][np.invert(f_drown.data[mask_eurat].mask)]))

#3. Natl
## DEF sea sst by domain
print('natl pt drown:',len(f_drown.data.data[mask_natl][np.invert(f_drown.data[mask_natl].mask)]))

#4. Hn
## DEF sea sst by domain
print('hn pt drown:',len(f_drown.data.data[mask_hn][np.invert(f_drown.data[mask_hn].mask)]))


#5. Hs
## DEF sea sst by domain
print('hs pt drown:',len(f_drown.data.data[mask_hs][np.invert(f_drown.data[mask_hs].mask)]))

#6. trop
## DEF sea sst by domain
print('trop pt drown:',len(f_drown.data.data[mask_trop][np.invert(f_drown.data[mask_trop].mask)]))

#7. glob
## DEF sea sst by domain
print('glob pt drown:',len(f_drown.data.data[mask_glob][np.invert(f_drown.data[mask_glob].mask)]))




## Read by computing diff  f_cmo - f & f_cmo - f_drown (Glob) & (lat>-60S)

# GLOB

# diff_glo = r_cmo.readfield(name[0]) - r.readfield(name_mer[0])
# len(diff_glo.data[np.invert(diff_glo.data.mask)]) # Retrouve nb c933 def


# diff_glo_taille=len(diff_glo.data[np.invert(diff_glo.data.mask)])
# print('diff_glo taille',diff_glo_taille)

# diff_glo_drown = f_cmo - f_drown
# diff_glo_drown_taille=len(diff_glo_drown.data[np.invert(diff_glo_drown.data.mask)])
# print('diff_glo taille',diff_glo_drown_taille)

# diff_glo_nord=diff_glo.data[mask_north]

# # =============================================================================
# # ##  Cartes sans forcer les valeurs def et non def 
# #     (fonctionnne car base sur masked epygram)
# # =============================================================================
# crs = ccrs.PlateCarree()
# # crs=ccrs.SouthPolarStereo(central_longitude=0.0, globe=None)
# fig, ax = f_diff.cartoplot(
# projection=crs, plot_method='scatter', colormap='seismic')
# #ax.set_extent([60, 85, 0, 30], ccrs.PlateCarree())
# ax.title.set_text(name[0]+'_'+str(expe)+
#   '_CMO_sansdrown_'+str(date)[0:10]+'-'+str(date)[11:13])
# #ax.set_extent([-110, -50, 40, 70], ccrs.PlateCarree())
# ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())




# ## Compare f et f drown
# f_c933=f-273.15
# crs = ccrs.PlateCarree()
# # crs=ccrs.SouthPolarStereo(central_longitude=0.0, globe=None)
# fig, ax = f_c933.cartoplot(
# projection=crs, plot_method='scatter', colormap='seismic')
# #ax.set_extent([60, 85, 0, 30], ccrs.PlateCarree())
# ax.title.set_text(name[0]+'_'+str(expe)+
#   '_x933_sansdrown_'+str(date)[0:10]+'-'+str(date)[11:13])
# #ax.set_extent([-110, -50, 40, 70], ccrs.PlateCarree())
# ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())


# f_c933_drown=f_drown-273.15
# crs = ccrs.PlateCarree()
# # crs=ccrs.SouthPolarStereo(central_longitude=0.0, globe=None)
# fig, ax = f_c933.cartoplot(
# projection=crs, plot_method='scatter', colormap='seismic')
# #ax.set_extent([60, 85, 0, 30], ccrs.PlateCarree())
# ax.title.set_text(name[0]+'_'+str(expe)+
#   '_c933_avecdrown_GMPE_'+str(date)[0:10]+'-'+str(date)[11:13])
# #ax.set_extent([-110, -50, 40, 70], ccrs.PlateCarree())
# ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())



### Cartes undef sans drown
# f_diff.data[mask_undef]=500
# f_diff.data[mask_def]=-999
# crs = ccrs.PlateCarree()
# # crs=ccrs.SouthPolarStereo(central_longitude=0.0, globe=None)
# fig, ax = f_diff.cartoplot(
# projection=crs, plot_method='scatter', colormap='seismic')
# #ax.set_extent([60, 85, 0, 30], ccrs.PlateCarree())
# ax.title.set_text(name[0]+'_'+str(expe)+
#   '_CMO_sansdrown_'+str(date)[0:10]+'-'+str(date)[11:13])
# #ax.set_extent([-110, -50, 40, 70], ccrs.PlateCarree())
# ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
# fig.savefig(chemin_undef+name[0]+'_'+str(expe)+'-'+'_'+str(date)
#     [0:10]+'-'+str(date)[11:13]+'CMO_mask_sans_drowning-'+cblock+'-'+ckind+str(term)+'.png')


# ### Cartes undef avec drown
# f_diff.data[mask_undef_drown]=500
# f_diff.data[mask_def_drown]=-999
# crs = ccrs.PlateCarree()
# # crs=ccrs.SouthPolarStereo(central_longitude=0.0, globe=None)
# fig, ax = f_diff.cartoplot(
# projection=crs, plot_method='scatter', colormap='seismic')
# #ax.set_extent([60, 85, 0, 30], ccrs.PlateCarree())
# ax.title.set_text(name[0]+'_'+str(expe) +'-'
#   '_Mask_avecdrown_GMPE_'+str(date)[0:10]+'-'+str(date)[11:13])
# #ax.set_extent([-110, -50, 40, 70], ccrs.PlateCarree())
# ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
# fig.savefig(chemin_undef+name[0]+'_'+str(expe)+'-'+'_'+str(date)
#     [0:10]+'-'+str(date)[11:13]+'-CMO_mask_avec_drowning-'+cblock+'-'+ckind+str(term)+'.png')

