#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 09:49:25 2023

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
import pandas as pd
import os
import numpy.ma as ma
epygram.init_env()

name_atm=[{'shortName':'slhf'}] # Chal Lat
#name_atm=[{'shortName':'sshf'}]
var='slhf'
#var='sshf'
date=datetime(2022,9,24,0)
MINMAX=[-200,200] #Sens
#MINMAX=[-400,400] #Lat

expe='GN84'
xpref='GKCJ'
cterm_for='96'
cterm_for_before='90'

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

def get_file_atm(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,'shouldfly':True,'cutoff':chain,'vapp':'arpege','vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr', 'origin':'historic', 'nativefmt':'grib','now':True} # nativefmt grib or fa
    fcref = toolbox.input(**fp)
    r=fcref[0].contents.data
    print('rrrr',r)
    return r

r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
r_for_atm_before=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for_before,geometry_atm,formats_atm,fill_atm)

r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
r_for_atm_ref_before=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for_before,geometry_atm,formats_atm,fill_atm)


x_expe=r_for_atm.readfield(name_atm[0])-r_for_atm_before.readfield(name_atm[0])
data=x_expe.getdata()
n,m=data.shape
tabval=np.zeros((n,m))
tabdiv=np.zeros((n,m))
tabval[:,:]=data
for k,j in np.ndindex(n,m):
   tabdiv[k,j]=np.divide(tabval[k,j],(6*3600))
x_expe.setdata(abs(tabdiv))



x_ref=r_for_atm_ref.readfield(name_atm[0])-r_for_atm_ref_before.readfield(name_atm[0])
data=x_ref.getdata()
n,m=data.shape
tabval=np.zeros((n,m))
tabdiv=np.zeros((n,m))
tabval[:,:]=data
for k,j in np.ndindex(n,m):
   tabdiv[k,j]=np.divide(tabval[k,j],(6*3600))
x_ref.setdata(abs(tabdiv))

crs=ccrs.PlateCarree()
#x_diff_eval.setdata(abs((r.readfield(name[0])-(r_mer.readfield(name_mer[0])-273.15)).getdata()))   ex abs
fig, ax = (x_expe-x_ref).cartoplot(projection=crs,plot_method='scatter',colormap='biais_lat_diff')

ax.set_extent([-90, -55, 10, 40], ccrs.PlateCarree())
#ax.title.set_text(str(name_atm)+'_Compar_'+str(expe)+'-'+str(xpref)+'_'+str(date)[0:10]+'-'+str(date)[11:13]+str(cterm_for))
plt.suptitle(str(name_atm)+'_Compar_abs_'+str(expe)+'_'+str(xpref)+'_'+str(date)[0:10]+'-'+str(date)[11:13]+str(cterm_for),fontsize=26)
plt.show()

fig.savefig('Compar_abs_'+str(expe)+'-'+str(xpref)+'_'+str(date)[0:10]+'-'+str(date)[11:13]+str(cterm_for)+'difftest_'+str(var)+'.watt.png')

