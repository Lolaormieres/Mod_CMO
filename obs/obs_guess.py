#!/usr/bin/env python
# -*- coding: utf-8 -*-

import usevortex
import os,sys
import shutil
import copy
from netCDF4 import Dataset
from datetime import datetime, date, timedelta
import numpy as np
import numpy.ma as ma
import scipy
import common
#import olive
import vortex
from vortex import tools
from vortex import toolbox
from footprints import proxy as fpx
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon
import matplotlib.path as mpath
import collections
import struct
#import epygram
#epygram.init_env()

###USEFULL OBJECTS
t = vortex.ticket()
e = t.env
sh = t.sh

begdat = tools.date.Date(2022,9,2,0)
enddat = tools.date.Date(2023,2,30,0)#31,21)#6,30,18)
dateList=[]
d=begdat

while (d <= enddat):
    dateList.append(d)
    d = d + tools.date.Period(hours=6)


xpid_dble = 'GN84'
xpid_ref = 'GKCJ'

kind = 'analyse'#'historic'#
block_can = 'canari'#'forecast'#
term = 6
geometry = vortex.data.geometries.get(tag='franmgsp')
fmt = 'fa'
co = 'assim'#'prod'#
vapp = 'arpege'
vconf = '4dvarfr'
model = 'surfex'#'arome'
filling = 'surf'#'atm'#
namespace = 'vortex.multi.fr'

list_fields = ['CLSTEMPERATURE','CLSHUMI.RELATIVE']
list_fields_sfx = []

obsvalue_tot = {}
varno_tot = {}
fg_depar_tot = {}
an_depar_tot = {}
anflag_tot = {}


for xp in [xpid_dble,xpid_ref]:
    obsvalue_tot[xp] = []
    varno_tot[xp] = []
    fg_depar_tot[xp] = []
    an_depar_tot[xp] = []
    anflag_tot[xp] = []

doGet = True


#baserequest_can = /home/martinezs/public/bin/odbsql "odbsql -q 'select lon, lat, date, time, statid, varno, obsvalue, fg_depar, an_depar, datum_anflag.final, from hdr,body, where varno=39 or varno=58 or varno=11 or varno=92' "
baserequest_can = "/home/martinezs/public/bin/odbsql -q 'select lon, lat, date, time, statid, varno, obsvalue, fg_depar, an_depar, datum_anflag.final, from hdr,body, where ((varno=39 or varno=58 or varno=11 or varno=92) and obstype=4)' "
# extraire obstype ... where obstype=4 (dribu)
#obstype=1 ou 4

nbobs = 0
for dat in dateList:

    #recupere les fichiers fa
    for xp in [xpid_dble,xpid_ref]:
        fpCan = {'experiment':xp,
              'kind':'observations',#'analysis',#'historic',
              'block':'surfan',#'forecast',#
              'role':'Observations',
              'format':'odb',
              'layout':'ecma',
              'part':'surf.tgz',
              'stage':'canari',
              'date':dat,#-epygram.util.datetime.timedelta(days=0, hours=6, minutes=0, seconds=0),
              'geometry':geometry,
              'local':'ECMA.[part]',
              'cutoff':co,
              'vapp':vapp,
              'vconf':vconf,
              'model':model,
              'namespace':namespace
              }

        fcCan = toolbox.input(**fpCan)[0]#usevortex.get_resources(getmode='epygram',**baseFootprintSurfanalyse)[0]#

        if (doGet):
            fcCan.get()

        oldpwd = os.getcwd()
        #os.system('rm -rf ECMA*')
        
        file_tar_can = (fcCan.locate()).split(';')[0]
        fcCan.locate()
        dir_file_can = file_tar_can.split('odb')[0]
        file_netcdf_can = file_tar_can[:-4]+".nc"
        os.chdir(dir_file_can)
        print (dir_file_can)
        #print (file_tar_can)
        os.system('ls -l')
        os.system('pwd')
        os.system('rm -rf ECMA*')
        

        ## Datarage 
        os.system('tar -xf '+file_tar_can)
        #os.system('ls')
        os.system('mv HIDDEN* ECMA_CAN')
        os.system(baserequest_can+" -i ECMA_CAN -f netcdf -o "+file_netcdf_can)
        os.system('rm -rf ECMA*')
        os.chdir(oldpwd)

        ncfile_can = Dataset(file_netcdf_can,'r')
        #os.system('ncdump -h '+file_netcdf_can)
        lon_can = ncfile_can.variables['col_1'][:]
        lat_can = ncfile_can.variables['col_2'][:]
        statid_can = ncfile_can.variables['col_5'][:]
        varno_can = ncfile_can.variables['col_6'][:]
        obsvalue_can = ncfile_can.variables['col_7'][:]
        fg_depar_can = ncfile_can.variables['col_8'][:]
        an_depar_can = ncfile_can.variables['col_9'][:]
        anflag_final_can = ncfile_can.variables['col_10'][:]
        ncfile_can.close()

        if (dat == begdat):
            obsvalue_tot[xp] = obsvalue_can
            varno_tot[xp] = varno_can
            fg_depar_tot[xp] = fg_depar_can
            an_depar_tot[xp] = an_depar_can
            anflag_tot[xp] = anflag_final_can
        else:
            obsvalue_tot[xp] = np.append(obsvalue_tot[xp],obsvalue_can)
            varno_tot[xp] = np.append(varno_tot[xp],varno_can)
            fg_depar_tot[xp] = np.append(fg_depar_tot[xp],fg_depar_can)
            an_depar_tot[xp] = np.append(an_depar_tot[xp],an_depar_can)
            anflag_tot[xp] = np.append(anflag_tot[xp],anflag_final_can)
        
    nbobs += np.shape(obsvalue_tot[xpid_dble])[0]




nb_obs_assim = {}
fg_min = {}
fg_max = {}
fg_mean = {}
fg_std = {}
an_min = {}
an_max = {}
an_mean = {}
an_std = {}

for xp in [xpid_dble,xpid_ref]:
    if (xp not in nb_obs_assim.keys()): nb_obs_assim[xp] = {}
    if (xp not in fg_min.keys()): fg_min[xp] = {}
    if (xp not in fg_max.keys()): fg_max[xp] = {}
    if (xp not in fg_mean.keys()): fg_mean[xp] = {}
    if (xp not in fg_std.keys()): fg_std[xp] = {}
    if (xp not in an_min.keys()): an_min[xp] = {}
    if (xp not in an_max.keys()): an_max[xp] = {}
    if (xp not in an_mean.keys()): an_mean[xp] = {}
    if (xp not in an_std.keys()): an_std[xp] = {}

    for varno in [39.,58.,11.,92.]:
        
        zz = np.where((varno_tot[xp] == varno) & (anflag_tot[xp] == 1.))[0]
        print (np.where(varno_tot[xp]))
        print (np.where(anflag_tot[xp] == 1.))
        nb_obs_assim[xp][varno] = np.shape(zz)[0]
        print (np.shape(zz))
       
        fg_min[xp][varno] = fg_depar_tot[xp][zz].min()
        fg_max[xp][varno] = fg_depar_tot[xp][zz].max()
        fg_mean[xp][varno] = fg_depar_tot[xp][zz].mean()
        fg_std[xp][varno] = fg_depar_tot[xp][zz].std()

        an_min[xp][varno] = an_depar_tot[xp][zz].min()
        an_max[xp][varno] = an_depar_tot[xp][zz].max()
        an_mean[xp][varno] = an_depar_tot[xp][zz].mean()
        an_std[xp][varno] = an_depar_tot[xp][zz].std()




namvar = {39.:'T2m',58.:'Hu2m',11.:'SST',92.:'Snow depth'}

for varno in [39.,58.,11.,92.]:
    print (namvar[varno],'fg departure')
    print ('nb = %6d (ref = %6d)' %(nb_obs_assim[xpid_dble][varno], nb_obs_assim[xpid_ref][varno]))
    print ('mean = %6.5f (%6.5f) std = %4.3f (%4.3f)' %(fg_mean[xpid_dble][varno], fg_mean[xpid_ref][varno], fg_std[xpid_dble][varno], fg_std[xpid_ref][varno]))
    print ('min = %6.5f (%6.5f) max = %4.3f (%4.3f)' %(fg_min[xpid_dble][varno], fg_min[xpid_ref][varno], fg_max[xpid_dble][varno], fg_max[xpid_ref][varno]))
    print (namvar[varno],'an departure')
    print ('nb = %6d (ref = %6d)' %(nb_obs_assim[xpid_dble][varno], nb_obs_assim[xpid_ref][varno]))
    print ('mean = %6.5f (%6.5f) std = %4.3f (%4.3f)' %(an_mean[xpid_dble][varno], an_mean[xpid_ref][varno], an_std[xpid_dble][varno], an_std[xpid_ref][varno]))
    print ('min = %6.5f (%6.5f) max = %4.3f (%4.3f)' %(an_min[xpid_dble][varno], an_min[xpid_ref][varno], an_max[xpid_dble][varno], an_max[xpid_ref][varno]))

        #mean_fg = fg_depar[np.where((anflag==1) | (anflag==4))[0]].mean()
        #mean_an = an_depar[np.where((anflag==1) | (anflag==4))[0]].mean()
        #std_fg = fg_depar[np.where((anflag==1) | (anflag==4))[0]].std()
        #std_an = an_depar[np.where((anflag==1) | (anflag==4))[0]].std()

        #mean_fg_used = fg_depar[np.where(anflag==1)[0]].mean()
        #mean_an_used = an_depar[np.where(anflag==1)[0]].mean()
        #std_fg_used = fg_depar[np.where(anflag==1)[0]].std()
        #std_an_used = an_depar[np.where(anflag==1)[0]].std()
        
        #fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2,ncols=2,figsize=(20,20))
        ##fig, axs = plt.subplots(2,2,sharey=True,sharex=True,tight_layout=True)
        #ax0.hist(fg_depar_ref[np.where(anflag_ref==1)[0]],bins=40,range=(-20.,20.))#,label=
        #ax0.set_title('REF fg_depar, mean='+str(mean_fg_ref_used)+', std='+str(std_fg_ref_used))
        #ax1.hist(an_depar_ref[np.where(anflag_ref==1)[0]],bins=40,range=(-20.,20.))#,label=
        #ax1.set_title('REF an_depar, mean='+str(mean_an_ref_used)+', std='+str(std_an_ref_used))
        #ax2.hist(fg_depar_sst[np.where(anflag_sst==1)[0]],bins=40,range=(-20.,20.))#,label=
        #ax2.set_title('SST fg_depar, mean='+str(mean_fg_sst_used)+', std='+str(std_fg_sst_used))
        #ax3.hist(an_depar_sst[np.where(anflag_sst==1)[0]],bins=40,range=(-20.,20.))#,label=
        #ax3.set_title('SST fg_depar, mean='+str(mean_an_sst_used)+', std='+str(std_an_sst_used))
        #fig.suptitle('SST histo nobs_ref='+str(np.shape(np.where(anflag_ref==1))[1])+', nobs_sst='+str(np.shape(np.where(anflag_sst==1))[1]))
        #fig.savefig('/cnrm/obs/data2/birmanc/surftemp/plots/histo_sst_used_'+xpid_sst+'_'+xpid_ref+'_'+str(dat)+'.png')

        #fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2,ncols=2,figsize=(20,20))
        ##fig, axs = plt.subplots(2,2,sharey=True,sharex=True,tight_layout=True)
        #ax0.hist(fg_depar_ref[np.where((anflag_ref==1) | (anflag_ref==4))[0]],bins=40,range=(-20.,20.))#,label=
        #ax0.set_title('REF fg_depar, mean='+str(mean_fg_ref)+', std='+str(std_fg_ref))
        #ax1.hist(an_depar_ref[np.where((anflag_ref==1) | (anflag_ref==4))[0]],bins=40,range=(-20.,20.))#,label=
        #ax1.set_title('REF an_depar, mean='+str(mean_an_ref)+', std='+str(std_an_ref))
        #ax2.hist(fg_depar_sst[np.where((anflag_sst==1) | (anflag_sst==4))[0]],bins=40,range=(-20.,20.))#,label=
        #ax2.set_title('SST fg_depar, mean='+str(mean_fg_sst)+', std='+str(std_fg_sst))
        #ax3.hist(an_depar_sst[np.where((anflag_sst==1) | (anflag_sst==4))[0]],bins=40,range=(-20.,20.))#,label=
        #ax3.set_title('SST fg_depar, mean='+str(mean_an_sst)+', std='+str(std_an_sst))
        #fig.suptitle('SST histo nobs_ref='+str(np.shape(np.where((anflag_ref==1) | (anflag_ref==4)))[1])+', nobs_sst='+str(np.shape(np.where((anflag_sst==1) | (anflag_sst==4)))[1]))
        #fig.savefig('/cnrm/obs/data2/birmanc/surftemp/plots/histo_sst_all_'+xpid_sst+'_'+xpid_ref+'_'+str(dat)+'.png')

        #fg_depar_tot[xp].extend(fg_depar[np.where((anflag==1) | (anflag==4))[0]])
        #an_depar_tot[xp].extend(an_depar[np.where((anflag==1) | (anflag==4))[0]])
        #anflag_tot[xp].extend(anflag[np.where((anflag==1) | (anflag==4))[0]])



'''
fg_depar_tot[xp] = np.array(fg_depar_tot[xp])
an_depar_tot[xp] = np.array(an_depar_tot[xp])
anflag_tot[xp] = np.array(anflag_tot[xp])


fig, ((ax0, ax1), (ax2, ax3), (ax4, ax5), (ax6, ax7)) = plt.subplots(nrows=4,ncols=2,figsize=(20,40))
#fig, axs = plt.subplots(2,2,sharey=True,sharex=True,tight_layout=True)
print(anflag_tot[xpid_ref],np.shape(anflag_tot[xpid_ref]))
print(np.where(anflag_tot[xpid_ref]==1.))
stop
print(fg_depar_tot[xpid_ref])
ax0.hist(fg_depar_tot[xpid_ref][np.where(anflag_tot[xpid_ref]==1)[0]],bins=40,range=(-20.,20.))#,label=
ax0.set_title('REF fg_depar, mean='+str(fg_depar_tot[xpid_ref][np.where(anflag_tot[xpid_ref]==1)[0]].mean())+', std='+str(fg_depar_tot[xpid_ref][np.where(anflag_tot[xpid_ref]==1)[0]].std()))
ax1.hist(an_depar_tot[xpid_ref][np.where(anflag_tot[xpid_ref]==1)[0]],bins=40,range=(-20.,20.))#,label=
ax1.set_title('REF an_depar, mean='+str(an_depar_tot[xpid_ref][np.where(anflag_tot[xpid_ref]==1)[0]].mean())+', std='+str(an_depar_tot[xpid_ref][np.where(anflag_tot[xpid_ref]==1)[0]].std()))

ax2.hist(fg_depar_tot[xpid_merc01][np.where(anflag_tot[xpid_merc01]==1)[0]],bins=40,range=(-20.,20.))#,label=
ax2.set_title('MERC1 fg_depar, mean='+str(fg_depar_tot[xpid_merc01][np.where(anflag_tot[xpid_merc01]==1)[0]].mean())+', std='+str(fg_depar_tot[xpid_merc01][np.where(anflag_tot[xpid_merc01]==1)[0]].std()))
ax3.hist(an_depar_tot[xpid_merc01][np.where(anflag_tot[xpid_merc01]==1)[0]],bins=40,range=(-20.,20.))#,label=
ax3.set_title('MERC1 an_depar, mean='+str(an_depar_tot[xpid_merc01][np.where(anflag_tot[xpid_merc01]==1)[0]].mean())+', std='+str(an_depar_tot[xpid_merc01][np.where(anflag_tot[xpid_merc01]==1)[0]].std()))

ax4.hist(fg_depar_tot[xpid_merc02][np.where(anflag_tot[xpid_merc02]==1)[0]],bins=40,range=(-20.,20.))#,label=
ax4.set_title('MERC2 fg_depar, mean='+str(fg_depar_tot[xpid_merc02][np.where(anflag_tot[xpid_merc02]==1)[0]].mean())+', std='+str(fg_depar_tot[xpid_merc02][np.where(anflag_tot[xpid_merc02]==1)[0]].std()))
ax5.hist(an_depar_tot[xpid_merc02][np.where(anflag_tot[xpid_merc02]==1)[0]],bins=40,range=(-20.,20.))#,label=
ax5.set_title('MERC2 an_depar, mean='+str(an_depar_tot[xpid_merc02][np.where(anflag_tot[xpid_merc02]==1)[0]].mean())+', std='+str(an_depar_tot[xpid_merc02][np.where(anflag_tot[xpid_merc02]==1)[0]].std()))

ax6.hist(fg_depar_tot[xpid_cep][np.where(anflag_tot[xpid_cep]==1)[0]],bins=40,range=(-20.,20.))#,label=
ax6.set_title('MERC2 fg_depar, mean='+str(fg_depar_tot[xpid_cep][np.where(anflag_tot[xpid_cep]==1)[0]].mean())+', std='+str(fg_depar_tot[xpid_cep][np.where(anflag_tot[xpid_cep]==1)[0]].std()))
ax7.hist(an_depar_tot[xpid_cer][np.where(anflag_tot[xpid_cep]==1)[0]],bins=40,range=(-20.,20.))#,label=
ax7.set_title('MERC2 an_depar, mean='+str(an_depar_tot[xpid_cep][np.where(anflag_tot[xpid_cep]==1)[0]].mean())+', std='+str(an_depar_tot[xpid_cep][np.where(anflag_tot[xpid_cep]==1)[0]].std()))

fig.suptitle('SST histo nobs_ref='+str(np.shape(np.where(anflag_tot[xpid_ref]==1))[1])+', nobs_merc1='+str(np.shape(np.where(anflag_tot[xpid_merc01]==1))[1])+', nobs_merc2='+str(np.shape(np.where(anflag_tot[xpid_merc02]==1))[1])+', nobs_cep='+str(np.shape(np.where(anflag_tot[xpid_cep]==1))[1]))
fig.savefig('/cnrm/obs/data2/birmanc/SST_mercator/plots/histo_sst_used_'+xpid_ref+'_'+xpid_merc01+'_'+xpid_merc02+'_'+xpid_cep+'_'+ref+str(begdat)+'_'+str(begdat)+'.png')
'''

