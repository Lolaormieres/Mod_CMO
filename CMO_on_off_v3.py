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
pathres = "/home/labadie/SCORES_CMO/scorestot/"
# experience à evaluer
expe = ['GL55', 'GN84','GO4A','GNYG','GKCJ']#OPER??
expe_mer = ['GL55', 'GN84','GO4A','GNYG','GN84']

# expe_mer=['GO48']
compar = 'prev-mer'+expe_mer[0]
# directory
Dir = "DF_SST_vf"
# Check if the directory exists
if not os.path.exists(Dir):
    # If it doesn't exist, create it
    os.makedirs(Dir)


# échéances
ech_fin = 102
ech_ini = 0
ech_step = 6

# nombre de niveaux
nb_level = 1

# MINMAX=[-2.5,2.5]
MINMAX=[-1.7,1.7]
MINMAX_sst=[18,33]
MINMAX=[-3e8,3e8]
# MINMAX=[-0.1e8,-3e8]
#MINMAX=[295,305]
# MINMAX_values=[700,1100]

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


# ATMO
name_atm=[{'shortName':'ssrd'},{'shortName':'10u'},{'shortName':'10v'}] #,{'shortNameECMF':'100v'}] 
name_atm=[{'shortName':'10u'}] 
name_atm=[{'shortName':'10v'}] 
#t temp level 0 level 0
# name_atm=[{'shortName':'t'}] 
# name: 'Relative humidity',
# shortName: 'r',
# name_atm=[{'shortName':'r'}] 
# name: '2 metre temperature',
# shortName: '2t',
# name_atm=[{'shortNameECMF':'2t'}] 
# name: '2 metre relative humidity',
# shortName: '2r',
# name_atm=[{'shortName':'2r'}] 
# name_atm=[{'shortName':'lcc'}] 
# short_name=['lcc'] 
# name_atm=[{'shortName':'prmsl'}] 
# name_atm=[{'shortName':'slhf'}] 
# short_name=['pres']
# short_name=['sst']

#ATMO FA
# name_atm=['SURFFLU.LAT.MEVA']
# name_atm=['SURFFLU.CHA.SENS']
# name_atm=['SURFFLU.LAT.MTOT']


#Diff journaliere
nb=1
date_init=datetime(2022,9,26,0)
date=date_init
dt_mer=timedelta(days=4)
dt=timedelta(days=1)

date_valid=date+dt_mer


expe='GKCJ'
expe2='GM6E'
expe3='GKPH'
xpref='GKCJ'

xmer ='Mercator'
#xpref='GKOR'
expe_drown='GMPE'
#geometry='global798c22'
geometry='global1798'
#term='24'
#chemin='/d0/images/ormieresl/'
chemin='/home/ormieresl/Scripts/Cartes'
chemin_pdg='/home/ormieresl/'
dir_images='/home/ormieresl/Routines/Storm_septembre/'


## ARP
cfields=''
cblock='surfan'
ckind='analysis'
cut='assim'
formats='fa'
fill_cmo='surf'
cterm='0'
#cblock='forecast'
#ckind='historic'
# historic.surfex.tl1798-c22+0078:00.fa # fichier Previ -

cblock_for='forecast'
ckind_for='historic'
cmodel='surfex'
cterm='0'
cmodel_for='surfex'
cut_for='production'
formats_for='fa'
fill_for=''
cterm_for='96'
cterm_for_before='90'

## Mer
#cblock_oc='c933_iter1_10_testITM'
cblock_oc='c933_iter1_10_testinverse'
cblock_oc='c933'
ckind_oc='geofields'
cfields_oc='ocean'
#cfields_oc='sst'
formats_oc='fa'
cmodel_oc='arpege'
cut_oc='assim'
fill_oc=''
cterm_oc='0'
## ARP Atmo
# grid.arpege-forecast.glob01+0090:00.grib # Fichier atm qu'on cherche à lire
fill_atm=''
grille="glob025"
geometry_atm=grille
cblock_atm='forecast'
nbmb=1
formats_atm='grib'
ckind_atm='gridpoint'
cblock_atm='forecast'
cfields_atm='grid'
cmodel_atm='arpege'
cut_atm='production'

## ARP fa altitude
# grid.arpege-forecast.glob01+0090:00.grib # Fichier atm qu'on cherche à lire
# fill_atm=''
# geometry_atm=geometry
# cblock_atm='forecast'
# ckind_atm='historic'
# nbmb=1
# formats_atm='fa'
# cblock_atm='forecast'
# cfields_atm=''
# cmodel_atm='arpege'
# cut_atm='production'

#%%
def get_file(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,
          'local':'tmpmer.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpmer.fic',openmode='r')
    os.remove('tmpmer.fic')
    return r

def get_file_atm(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,
          'local':'tmpatm.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr', 'origin':'historic',
          'nativefmt':'grib','now':True} # nativefmt grib or fa
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpatm.fic',openmode='r')
    return r
   
#%%
## DIFF EXPE
def plot_diffmean():
    for i in range(len(name)):
        print(i,name[i])
        crs=ccrs.PlateCarree()
        #fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
        # fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
       
        fig, ax = x_diff.cartoplot(projection=crs,plot_method='scatter',colormap='seismic') #'biais_v4_grey'(len(name)):

        ax.set_extent([-95, -60, 10, 40], ccrs.PlateCarree())
        #ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
        # ax.set_extent([-85, -50, 10, 40], ccrs.PlateCarree())
        ax.title.set_text(str(name_atm)+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        # ax.title.set_text(name[i]+'_'+str(xpref)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        
        #fig.savefig(dir_images+name[i]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-'+cblock+'-'+ckind+str(cterm)+'.png')

def plot_values():
    for i in range(len(name)):
        crs=ccrs.PlateCarree()
        # fig, ax = x.cartoplot(projection=crs,plot_method='scatter',colormap='Blues') #'biais_v4_grey'(len(name)):
        fig, ax = x.cartoplot(projection=crs,plot_method='scatter',colormap='Blues') #'biais_v4_grey'(len(name)):
        ax.quiver()
        # ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
        ax.set_extent([-95, -60, 10, 40], ccrs.PlateCarree())
        ax.title.set_text(str(name_atm)+'_'+str(expe)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        # ax.title.set_text(short_name[0]+'_'+str(expe)+' '+str(date)[0:10]+'-'+str(date)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-P0')

        # ax.title.set_text(name[i]+'_'+str(xpref)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        
        #fig.savefig(dir_images+short_name[0]+'_'+str(expe)+'-'+'_'+str(date)[0:10]+'-'+str(date)[11:13]+'_'+str(cterm_for)+'_P0.png')
        

# =============================================================================
# #DEUX EXPE
# =============================================================================

# for t in range(nb):
#     print(date)
    
#     try:

#         # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
#         r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
#         r_for_ref=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
          
#         # r_mer = get_file(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
#         if t == 0:
#           x_diff = r_for.readfield(name_sst[0])-(r_for_ref.readfield(name_sst[0]))
          
#         else:
#           x_diff=r_for.readfield(name_sst[0])-(r_for_ref.readfield(name_sst[0]))+x_diff
#     except:
#         print('Pas de fichier')
#     date=date+dt
    
# x_diff.setdata(x_diff.getdata()/nb)
# plot_diffmean()


# =============================================================================
# #Value EXPE
# =============================================================================
# for t in range(nb):
#     print(date)
    
#     try:
#         # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
        
#         r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
#         # r_for_ref=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
   
#         #ATM:
#         #r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#         # r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#           # r_mer = get_file(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
#         # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
#         # r_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel, cut, cterm)
       
#         x = r_for.readfield(name_sst[0])
#         X = r_for.geometry.get_lonlat_grid()  # lon et lat
#         lon = x[0]
#         lat = X[1]
        
#         mask=(lon >-80) & (lon <-50) & (lat>10) & (lat>40)
#         # x = r_for_atm_ref.readfield(name_atm[0])
#     except:
#         print('Pas de fichier')

#     # x.setdata(x)
#     plot_values()
#     date=date+dt

    # print('miiiiiiiiiiin',np.nanmin(x.data[mask]))
    # print('maaaaaaaaaaax',np.nanmax(x.data[mask]))


# =============================================================================
# #MERCATOR  
# =============================================================================
# for t in range(nb):
#     print(date)
    
#     try:
#        # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
        # r_for=get_file_atm(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
        # r_mer = get_file_atm(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
#         if t == 0:
#             x_diff = r_for.readfield(name[0])-(r_mer.readfield(name_mer[0]))
#         else:
#             x_diff=r_for.readfield(name[0])-(r_mer.readfield(name_mer[0]))+x_diff
#     except:
#         print('Pas de fichier')
#     date=date+dt
    
# x_diff.setdata(x_diff.getdata()/nb)
# plot_diffmean()

# =============================================================================
# SST VALUES no function
# =============================================================================
print('********************************cterm_for=',cterm_for)
#r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
# r_mer = get_file(expe,date_valid,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
# r_mer_before = get_file(expe,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)


r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
#r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
# r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm) 
print('********************************cterm_for_b=',cterm_for_before)
# r_for_before=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
#x_diff = r_for_atm.readfield(name_atm[0])
  # - r_for_atm_ref.readfield(name_atm[0])
 
X=r_for.geometry.get_lonlat_grid()
lon_expe=X[0]
lat_expe=X[1]

x_diff=r_for.readfield(name_sst[0]) - 273.15
# x_u=r_for_atm.readfield(name_u[0]) 
# x_v=r_for_atm.readfield(name_v[0]) 
        #r_mer.readfield(name_mer[0])) - 273.15
# - r_for_before.readfield(name_sst[0])
# x_diff2=r_for_before

mask=(lon_expe<-55) & (lon_expe>-90) & (lat_expe>10) & (lat_expe<40)
Minimum=x_diff.data[mask].min()
Maximum=x_diff.data[mask].max()
Moy    =x_diff.data[mask].mean()

crs=ccrs.PlateCarree()
# fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
   
fig, ax = x_diff.cartoplot(projection=crs,plot_method='scatter',minmax=MINMAX_sst,colormap='biais_sst') #'biais_v4_grey', 'biais_sst_cyclone_thiner', 'biais_sst_v1'
# fig, ax = x_u.cartoplot(projection=crs,plot_method='scatter',minmax=MINMAX_sst,colormap='biais_sst') 

ax.set_extent([-90, -55, 10, 40], ccrs.PlateCarree())
#ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
# ax.title.set_text(str(name_sst)+'_'+str(expe)+'_Date_valid'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'P'+str(cterm_for))
plt.suptitle(str(name_sst)+'_'+str(expe)+'_Date_valid: '+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_P'+str(cterm_for),fontsize=26)
plt.show()

# #fig.savefig(dir_images+name[i]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-'+cblock+'-'+ckind+str(cterm)+'.png')


# # =============================================================================
# # SST diff no function
# # =============================================================================
# print('********************************cterm_for=',cterm_for)
# r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
# #r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
# print('********************************cterm_for_b=',cterm_for_before)
# r_for_before=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for_before,geometry,formats,fill_cmo)
# #x_diff = r_for_atm.readfield(name_atm[0])
#   # - r_for_atm_ref.readfield(name_atm[0])
 
 
# X=r_for.geometry.get_lonlat_grid()
# lon_expe=X[0]
# lat_expe=X[1]

# x_diff=r_for.readfield(name_sst[0]) - r_for_before.readfield(name_sst[0])
# x_diff2=r_for_before


# mask=(lon_expe<-65) & (lon_expe>-95) & (lat_expe>10) & (lat_expe<40)
# Minimum=x_diff.data[mask].min()
# Maximum=x_diff.data[mask].max()
# Moy    =x_diff.data[mask].mean()

# crs=ccrs.PlateCarree()
# #fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
# # fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
   
# fig, ax = x_diff.cartoplot(projection=crs,plot_method='scatter',colormap='seismic') #'biais_v4_grey'(len(name)), 'biais_sst_cyclone_thiner'

# ax.set_extent([-90, -55, 10, 40], ccrs.PlateCarree())

# #ax.title.set_text(str(name)+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
# # ax.title.set_text(name[i]+'_'+str(xpref)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
# plt.suptitle(str(name_sst)+'_'+str(expe)+'_Date_valid: '+str(date_valid)[0:10]+'-'+str(date_valid)[11:13]+'_P96-P0',fontsize=26)
# plt.show()
# # #fig.savefig(dir_images+name[i]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-'+cblock+'-'+ckind+str(cterm)+'.png')



# =============================================================================
# SANS FUNCTION ATM
# =============================================================================
# # print('********************************cterm_for=',cterm_for,'name_atm[0]=',name_atm[0])
# # r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
# # #r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
# # print('********************************cterm_for_b=',cterm_for_before)
# # r_for_atm_before=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for_before,geometry_atm,formats_atm,fill_atm)

# #x_diff = r_for_atm.readfield(name_atm[0])
#  # - r_for_atm_ref.readfield(name_atm[0])
 
# x_diff=r_for_atm
# x_diff2=r_for_atm_before

# crs=ccrs.PlateCarree()
# #fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
# # fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
   
# fig, ax = x_diff.readfield(name_atm[0]).cartoplot(projection=crs,plot_method='scatter',colormap='seismic') #'biais_v4_grey'(len(name)):

# ax.set_extent([-85, -50, 10, 40], ccrs.PlateCarree())
# #ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
# # ax.set_extent([-85, -50, 10, 40], ccrs.PlateCarree())
# ax.title.set_text(str(name_atm)+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
# # ax.title.set_text(name[i]+'_'+str(xpref)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')

# #fig.savefig(dir_images+name[i]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-'+cblock+'-'+ckind+str(cterm)+'.png')


# =============================================================================
# # ATMOSPHERE
# =============================================================================

# # Diff deux exp
# for t in range(nb):
#     print(date)
    
#     try:
#         # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
#         # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
#         # r_for_ref=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
        
#         #ATM:
#         r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#         # r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#         r_for_atm_before=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for_before,geometry_atm,formats_atm,fill_atm)
        
#         # r_mer = get_file(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
#         # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
#         # r_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel, cut, cterm)
#         if t == 0:
#             x_diff = r_for_atm.readfield(name_atm[0])-(r_for_atm_before.readfield(name_atm[0]))
#         else:
#             x_diff=r_for_atm.readfield(name_atm[0])-(r_for_atm_before.readfield(name_atm[0]))+x_diff
#     except:
#         print('Pas de fichier')
#     date=date+dt
    
# x_diff.setdata(x_diff.getdata()/nb)
# plot_diffmean()


# =============================================================================
# Values ATM
# =============================================================================


# for t in range(nb):
#     print(date)
    
#     try:
#         # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
        
#         # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
#         # r_for_ref=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
        
#         #ATM:
#         r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#         r_for_atm_before=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for_before,geometry_atm,formats_atm,fill_atm)
       
#         # r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#           # r_mer = get_file(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
#         # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
#         # r_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel, cut, cterm)
#         # x = r_for_atm.readfield(name_atm[0]) 
#         x = (r_for_atm.readfield(name_atm[0]))
#           # - r_for_atm_ref.readfield(name_atm[0])
 
#         X=r_for_atm.geometry.get_lonlat_grid()
#         lon_expe=X[0]
#         lat_expe=X[1]
        
        
#         mask=(lon_expe<-65) & (lon_expe>-95) & (lat_expe>10) & (lat_expe<40)
#         Minimum=x.data[mask].min()
#         print('min',Minimum)
#         Maximum=x.data[mask].max()
#         print('max',Maximum)
#         Moy    =x.data[mask].mean()
#         print('moy', Moy)
        
#               #- r_for_atm_before.readfield(name_atm[0]) ) 
#         # x2 = x.data.data/(6*3600)
#         # x = r_for_atm_ref.readfield(name_atm[0])
#     except:
#         print('Pas de fichier')
#     date=date+dt
    
#     # x.setdata(x)
#     plot_values()


# =============================================================================
# # VALUES ATM no func
# 
# =============================================================================
# for t in range(nb):
#     print(date)
    

# r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)

# r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
# r_for_ref=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)

# #ATM:
# r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
# r_for_atm_before=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for_before,geometry_atm,formats_atm,fill_atm)
   
# # r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#   # r_mer = get_file(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
# # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
# # r_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel, cut, cterm)
# # x = r_for_atm.readfield(name_atm[0]) 
# x = (r_for_atm.readfield(name_atm[0]))
#   # - r_for_atm_ref.readfield(name_atm[0])
 
# # X=r_for_atm.geometry.get_lonlat_grid()
# # lon_expe=X[0]
# # lat_expe=X[1]


# # mask=(lon_expe<-65) & (lon_expe>-95) & (lat_expe>10) & (lat_expe<40)
# # Minimum=x.data[mask].min()
# # print('min',Minimum)
# # Maximum=x.data[mask].max()
# # print('max',Maximum)
# # Moy    =x.data[mask].mean()
# # print('moy', Moy)

#       #- r_for_atm_before.readfield(name_atm[0]) ) 
# # x2 = x.data.data/(6*3600)
# # x = r_for_atm_ref.readfield(name_atm[0])

# # date=date+dt




# # fig
# crs=ccrs.PlateCarree()
# # fig, ax = x.cartoplot(projection=crs,plot_method='scatter',colormap='Blues') #'biais_v4_grey'(len(name)):
# fig, ax = x.cartoplot(projection=crs,plot_method='scatter',minmax=MINMAX,colormap='BrBG') #'biais_v4_grey'(len(name)):
# # ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
# ax.set_extent([-95, -60, 10, 40], ccrs.PlateCarree())
# # ax.title.set_text(str(name_atm)+'_'+str(expe)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for))
# # ax.title.set_text(short_name[0]+'_'+str(expe)+' '+str(date)[0:10]+'-'+str(date)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-P0')
# plt.suptitle(str(name_atm)+'_'+str(expe)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for),fontsize=26)
# plt.show()
# # ax.title.set_text(name[i]+'_'+str(xpref)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        
    