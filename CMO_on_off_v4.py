#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 17:33:37 2023

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
dir_im = '/cnrm/recyf/Data/users/ormieresl/plot_sst_cyclones/'

# experience à evaluer
expe = ['GL55', 'GN84', 'GO4A', 'GNYG', 'GKCJ']  # OPER??
expe_mer = ['GL55', 'GN84', 'GO4A', 'GNYG', 'GN84']

# expe_mer=['GO48']
compar = 'prev-mer'+expe_mer[0]
# directory
Dir = "Plot_cartes_cyclones"
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
MINMAX = [-1.7, 1.7]
MINMAX_sst = [18, 33]
MINMAX = [-3e8, 3e8]
# MINMAX=[-0.1e8,-3e8]
# MINMAX=[295,305]
# MINMAX_values=[700,1100]




geometry = 'global1798'



# lecture des variables
name_sst = ['SFX.SST']
name = ['SFX.TEMP_OC'+str(i) for i in range(1, nb_level+1)]
name_sal = ['SFX.SALT_OC'+str(i) for i in range(1, nb_level+1)]
name_mer = ['SURF.THETAO'+str(i) for i in range(1, nb_level+1)]
name_sal_mer = ['SURF.SALINO'+str(i) for i in range(1, nb_level+1)]
name_u = ['SFX.UCUR_OC'+str(i) for i in range(1, nb_level+1)]
name_v = ['SFX.VCUR_OC'+str(i) for i in range(1, nb_level+1)]

# ATMO
name_atm = [{'shortName': 'ssrd'}, {'shortName': '10u'},
            {'shortName': '10v'}]  # ,{'shortNameECMF':'100v'}]
name_atm = [{'shortName': '10u'}]
# name_atm=[{'shortName':'10v'}]
# t temp level 0 level 0
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
name_atm=[{'shortName':'prmsl'}]
name_atm_lat=[{'shortName':'slhf'}]
# short_name=['pres']
# short_name=['sst']

# ATMO FA
# name_atm=['SURFFLU.LAT.MEVA']
# name_atm=['SURFFLU.CHA.SENS']
# name_atm=['SURFFLU.LAT.MTOT']


# ARP
cfields = ''
cblock = 'surfan'
ckind = 'analysis'
cut = 'assim'
formats = 'fa'
fill_cmo = 'surf'
cterm = '0'
# cblock='forecast'
# ckind='historic'
# historic.surfex.tl1798-c22+0078:00.fa # fichier Previ -

cblock_for = 'forecast'
ckind_for = 'historic'
cmodel = 'surfex'
cterm = '0'
cmodel_for = 'surfex'
cut_for = 'production'
formats_for = 'fa'
fill_for = ''
cterm_for = '96'
cterm_for_before = '0'

# Mer
# cblock_oc='c933_iter1_10_testITM'
cblock_oc = 'c933_iter1_10_testinverse'
cblock_oc = 'c933'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
# cfields_oc='sst'
formats_oc = 'fa'
cmodel_oc = 'arpege'
cut_oc = 'assim'
fill_oc = ''
cterm_oc = '0'
# ARP Atmo
# grid.arpege-forecast.glob01+0090:00.grib # Fichier atm qu'on cherche à lire
fill_atm = ''
grille = "glob025"
geometry_atm = grille
cblock_atm = 'forecast'
nbmb = 1
formats_atm = 'grib'
ckind_atm = 'gridpoint'
cblock_atm = 'forecast'
cfields_atm = 'grid'
cmodel_atm = 'arpege'
cut_atm = 'production'


def get_file(xpid, date, rep, typfic, fields, mod, chain, ech, geo, formats, fill):
    # MODEL
    fp = {'experiment': xpid, 'kind': typfic, 'block': rep, 'fields': fields,
          'date': date, 'term': ech, 'geometry': geo, 'format': formats, 'filling': fill,
          'local': 'tmpcartes.fic', 'cutoff': chain, 'vapp': 'arpege',
          'vconf': '4dvarfr', 'model': mod, 'namespace': 'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpcartes.fic', openmode='r')
    os.remove('tmpcartes.fic')
    return r


def get_file_atm(xpid, date, rep, typfic, fields, mod, chain, ech, geo, formats, fill):
    # MODEL
    fp = {'experiment': xpid, 'kind': typfic, 'block': rep, 'fields': fields,
          'date': date, 'term': ech, 'geometry': geo, 'format': formats, 'filling': fill,
          'local': 'tmpatm.fic', 'cutoff': chain, 'vapp': 'arpege',
          'vconf': '4dvarfr', 'model': mod, 'namespace': 'vortex.multi.fr', 'origin': 'historic',
          'nativefmt': 'grib', 'now': True}  # nativefmt grib or fa
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpatm.fic', openmode='r')
    return r


# =============================================================================
# SST VALUES no function
# =============================================================================
# -------------------------------------------------------------------------------------------------
# date
# -------------------------------------------------------------------------------------------------
cterm = '0'
cterm_for = '96'
cterm_for_before = '0'

expe = 'GOJQ'
xpref = 'GOJQ'

nb = 1
date = datetime(2022, 9, 28, 0)
dt_valid = timedelta(hours=int(cterm_for_before))
# dt=timedelta(days=1)
date_valid = date+dt_valid

P=cterm_for

print('********************************cterm=', cterm_for)
# r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
#r_mer = get_file(expe,date_valid,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)

print('********************************cterm_for=', cterm_for)
r_for = get_file(expe, date, cblock_for, ckind_for, cfields,cmodel, cut_for, cterm_for, geometry, formats, fill_cmo)
r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
# r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)

print('********************************cterm_for_before=', cterm_for_before)
# r_for_before = get_file(expe, date, cblock_for, ckind_for, cfields,
#                         cmodel, cut_for, cterm_for_before, geometry, formats, fill_cmo)
# #r_mer_before = get_file(expe,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
# r_for_atm_before=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for_before,geometry_atm,formats_atm,fill_atm)


# ------------------------------------------------------------------------------
X = r_for.geometry.get_lonlat_grid()  # Not for ATM
lon_expe = X[0]
lat_expe = X[1]

mask= (lon_expe < -55) & (lon_expe > -90) & (lat_expe > 10) & (lat_expe < 40)


x = r_for.readfield(name_sst[0]) - 273.15
x_atm = r_for_atm.readfield(name_atm[0])


## PLOT SST
# ------------------------------------------------------------------------------
crs = ccrs.PlateCarree()
crs = ccrs.LambertConformal()
# fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
# levels=np.arange(860,1100,20)
fig, ax = x.cartoplot(projection=crs, plot_method='contourf', colormap='biais_sst')

# pas = 300
# pmin=np.min(x_atm.getdata())
# pmax=np.max(x_atm.getdata())
# levels=list(np.arange(101000+pas,np.min(x_atm.getdata()),-pas))
# levels.reverse()
# numref=len(levels)-1
# levels += list(np.arange(101000+2*pas,np.max(x_atm.getdata()), pas))

# if pmax > levels[-1]: levels += [levels[-1] + pas]
# lw=[2.]*len(levels)
# lw[numref]=3.5
# labs=dict(manual=False,inline=True,fmt=lambda x:' {:.0f} '.format(x/100))
# x_atm.cartoplot(ax=ax,plot_method='contour',contour_kw=dict(colors=['brown'],levels=levels, linewidths=lw),  fig = fig)

# ax.set_extent([-90, -55, 10, 40], ccrs.PlateCarree())
# ax.set_extent([-90, -55, 10, 40], ccrs.PlateCarree())
x_atm.cartoplot(ax=ax, projection=crs,plot_method='contour',  fig = fig, linewidth=0.5)
ax.set_extent([-90, -55, 10, 40], ccrs.PlateCarree())

ax.tick_params(axis='y', labelsize=15)
# # clb.ax.set_title('Bias (°C)',fontsize=19)
plt.yticks(fontsize=20)  


# # taken fromd undef_zoom tickbar --> make bigger lon/lat
# gl = ax.gridlines(ccrs.LambertConformal(), draw_labels=True,
#                 linewidth=2, color='gray', alpha=0.5, linestyle='--')
# # gl.xlabels_top = None
# # gl.ylabels_left = True
# # gl.xlabels_bottom = True
# # gl.ylabels_right = None
# gl.xlabel_style = {'size': 17, 'color': 'black', 'weight':'bold'}
# gl.ylabel_style = {'size': 17, 'color': 'black','weight':'bold'}

# ax.tick_params(axis='y', labelsize=18)
# plt.yticks(fontsize=20) 
# ax.tick_params(axis='x', labelsize=18)
# plt.show()

# fig, ax = x_atm.cartoplot(projection=crs,plot_method='contour')

# clabel_kw=labs
# Colorbar
#'biais_v4_grey', 'biais_sst_cyclone_thiner', 'biais_sst_v1','biais_sst_cyclone_large','biais_sst_cyclone_restricted'


#ax.set_extent([-90, -55, 10, 40], ccrs.PlateCarree())
plt.suptitle(str(name_sst[0])+'_'+str(expe)+'_Date_valid: '+str(date_valid)
             [0:10]+'-'+str(date_valid)[11:13]+'_P'+str(P), fontsize=26)
plt.show()

figname=dir_im+str(name_sst[0])+'_'+str(expe)+'_'+str(date)[0:10]+'-'+str(date)[11:13]+'_P'+str(P)+'.png'
fig.savefig(figname,dpi=200, format='png',bbox_inches=None, pad_inches=0.1)
figname=dir_im+str(name_sst[0])+'_'+str(expe)+'_'+str(date)[0:10]+'-'+str(date)[11:13]+'_P'+str(P)+'.pdf'
fig.savefig(figname,dpi=200, format='pdf',bbox_inches=None, pad_inches=0.1)





#x_atm.cartoplot(ax=ax,fig=fig,plot_method='contour')

# # #fig, ax = x_atm.cartoplot(projection=crs,plot_method='scatter')   
# # fig, ax = x_atm.cartoplot(
# #     projection=crs, plot_method='scatter')          
# # ax.set_extent([-90, -55, 10, 40], ccrs.PlateCarree())
# titre = str(name_sst[0])+'_'+str(expe)+'_Date_valid: '+str(date_valid)[0:10]+'-'+str(date_valid)[11:13]+'_P'+str(P)
# rcparams=[(('font',), {'family': 'serif'}),(('axes',), {'titlepad': 24.})]
# pas = 300
#     #champ = grib.readfield({'shortName':'prmsl'}).extract_zoom(zoom)
# pmin=np.min(x_atm.getdata())
# pmax=np.max(x_atm.getdata())
# levels=list(np.arange(101000+pas,np.min(x_atm.getdata()),-pas))
# levels.reverse()
# numref=len(levels)-1
# levels += list(np.arange(101000+2*pas,np.max(x_atm.getdata()), pas))
# if pmax > levels[-1]: levels += [levels[-1] + pas]
# lw=[2.]*len(levels)
# lw[numref]=3.5
# labs=dict(manual=False,inline=True,fmt=lambda x:' {:.0f} '.format(x/100))
# x_atm.cartoplot(ax=ax,plot_method='contour',clabel_kw=labs,contour_kw=dict(colors=['brown'],levels=levels, linewidths=lw), contourlabel=True, fig = fig, title = titre, rcparams = rcparams)
         
# ###
# sst_nc=x_diff.dump_to_nc('sst.nc',variablename='SFX.SST')
