#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:30:16 2023

@author: ormieresl
"""
#! /usr/bin/python

# coding: utf-8

#import bronx
from matplotlib import font_manager as fm, rcParams
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

#name = ['SFX.SST']
#name_mer = ['SURFSST.CLIM.']
name = ['SFX.TEMP_OC1']  # , 'SFX.TEMP_OC2', 'SFX.TEMP_OC3', 'SFX.TEMP_OC4', 'SFX.TEMP_OC5', 'SFX.TEMP_OC6', 'SFX.TEMP_OC7', 'SFX.TEMP_OC8', 'SFX.TEMP_OC9', 'SFX.TEMP_OC10', 'SFX.TEMP_OC11', 'SFX.TEMP_OC12']#,
# 'SFX.TEMP_OC13','SFX.TEMP_OC14', 'SFX.TEMP_OC15', 'SFX.TEMP_OC16', 'SFX.TEMP_OC17', 'SFX.TEMP_OC18', 'SFX.TEMP_OC19', 'SFX.TEMP_OC20', 'SFX.TEMP_OC21', 'SFX.TEMP_OC22', 'SFX.TEMP_OC23', 'SFX.TEMP_OC24', 'SFX.TEMP_OC25', 'SFX.TEMP_OC26']
name_sal = ['SFX.SALT_OC1', 'SFX.SALT_OC2', 'SFX.SALT_OC3', 'SFX.SALT_OC4', 'SFX.SALT_OC5', 'SFX.SALT_OC6', 'SFX.SALT_OC7', 'SFX.SALT_OC8', 'SFX.SALT_OC9', 'SFX.SALT_OC10', 'SFX.SALT_OC11', 'SFX.SALT_OC12', 'SFX.SALT_OC13',
            'SFX.SALT_OC14', 'SFX.SALT_OC15', 'SFX.SALT_OC16', 'SFX.SALT_OC17', 'SFX.SALT_OC18', 'SFX.SALT_OC19', 'SFX.SALT_OC20', 'SFX.SALT_OC21', 'SFX.SALT_OC22', 'SFX.SALT_OC23', 'SFX.SALT_OC24', 'SFX.SALT_OC25', 'SFX.SALT_OC26']
# , 'SURF.THETAO2', 'SURF.THETAO3', 'SURF.THETAO4', 'SURF.THETAO5', 'SURF.THETAO6', 'SURF.THETAO7', 'SURF.THETAO8', 'SURF.THETAO9', 'SURF.THETAO10', 'SURF.THETAO11', 'SURF.THETAO12']#,
name_mer = ['SURF.THETAO1']
# 'SURF.THETAO13','SURF.THETAO14', 'SURF.THETAO15', 'SURF.THETAO16', 'SURF.THETAO17', 'SURF.THETAO18', 'SURF.THETAO19', 'SURF.THETAO20', 'SURF.THETAO21', 'SURF.THETAO22', 'SURF.THETAO23', 'SURF.THETAO24', 'SURF.THETAO25', 'SURF.THETAO26']
name_sal_mer = ['SURF.SALINO1', 'SURF.SALINO2', 'SURF.SALINO3', 'SURF.SALINO4', 'SURF.SALINO5', 'SURF.SALINO6', 'SURF.SALINO7', 'SURF.SALINO8', 'SURF.SALINO9', 'SURF.SALINO10', 'SURF.SALINO11', 'SURF.SALINO12', 'SURF.SALINO13',
                'SURF.SALINO14', 'SURF.SALINO15', 'SURF.SALINO16', 'SURF.SALINO17', 'SURF.SALINO18', 'SURF.SALINO19', 'SURF.SALINO20', 'SURF.SALINO21', 'SURF.SALINO22', 'SURF.SALINO23', 'SURF.SALINO24', 'SURF.SALINO25', 'SURF.SALINO26']
name_u = ['SFX.UCUR_OC1', 'SFX.UCUR_OC2', 'SFX.UCUR_OC3', 'SFX.UCUR_OC4', 'SFX.UCUR_OC5', 'SFX.UCUR_OC6', 'SFX.UCUR_OC7', 'SFX.UCUR_OC8', 'SFX.UCUR_OC9', 'SFX.UCUR_OC10', 'SFX.UCUR_OC11', 'SFX.UCUR_OC12', 'SFX.UCUR_OC13',
          'SFX.UCUR_OC14', 'SFX.UCUR_OC15', 'SFX.UCUR_OC16', 'SFX.UCUR_OC17', 'SFX.UCUR_OC18', 'SFX.UCUR_OC19', 'SFX.UCUR_OC20', 'SFX.UCUR_OC21', 'SFX.UCUR_OC22', 'SFX.UCUR_OC23', 'SFX.UCUR_OC24', 'SFX.UCUR_OC25', 'SFX.UCUR_OC26']
name_v = ['SFX.VCUR_OC1', 'SFX.VCUR_OC2', 'SFX.VCUR_OC3', 'SFX.VCUR_OC4', 'SFX.VCUR_OC5', 'SFX.VCUR_OC6', 'SFX.VCUR_OC7', 'SFX.VCUR_OC8', 'SFX.VCUR_OC9', 'SFX.VCUR_OC10', 'SFX.VCUR_OC11', 'SFX.VCUR_OC12', 'SFX.VCUR_OC13',
          'SFX.VCUR_OC14', 'SFX.VCUR_OC15', 'SFX.VCUR_OC16', 'SFX.VCUR_OC17', 'SFX.VCUR_OC18', 'SFX.VCUR_OC19', 'SFX.VCUR_OC20', 'SFX.VCUR_OC21', 'SFX.VCUR_OC22', 'SFX.VCUR_OC23', 'SFX.VCUR_OC24', 'SFX.VCUR_OC25', 'SFX.VCUR_OC26']


# Lecture des niveaux de la CMO
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
ech = ''
# Mercator
cblock_oc = 'c933'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
#cfields_oc = 'sst'
cmodel_oc = 'arpege'
cut_oc = 'assim'

# geometry='global798c22'
geometry = 'global1798'
# chemin='/d0/images/ormieresl/'
chemin = '/home/ormieresl/Scripts/Dataframe_mer/'
ech = cterm
Dir = '/home/ormieresl/Routines/DF_bouees_ech_vf/'
# Dir_im='DF_bouees_ech_vf_new/plot_ech/'
Dir_im = '/cnrm/recyf/Data/users/ormieresl/plot_ech_bouees/'

# Definition color
color_GL55 = '#648FFF'  # GL55
color_GKOR = '#785EF0'  # GKOR
color_GKPH = 'red'  # GKPH '#DC267F'
color_GM6E = '#FE6100'  # GM6E
color_GMOT = '#009E73'  # GMOT
color_GK7C = 'dimgrey'
color_GN3C = 'gold'
color_GN84 = 'lightseagreen'
color_GNBA = 'midnightblue'


# Definition color Pale
color_GL55 = 'steelblue'  # GL55 '#648FFF'
color_GKOR = 'hotpink'  # GKOR
color_GKPH = 'mediumpurple'  # GKPH '#DC267F'
color_GM6E = 'paleturquoise'
color_GMOT = 'peachpuff'
color_GN84 = 'r'
color_GO4A = 'rosybrown'
color_GNSR = 'palegreen'
color_GNYG = 'yellow'

color_expe = color_GM6E
color_mer = 'gold'

# %%
# =============================================================================
# ## Lectur ficher fa
# =============================================================================


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
    # return fcref.contents.data
    # fcref.container.filename
    return r
# %%


zone = ['eurat', 'nordat', "tropiques", "hn20", "hs20", "glob", "med"]
expe = ['GL55', 'GKOR', 'GKPH', 'GN84', 'GO4A']
color = [color_GL55, color_GKOR, color_GKPH, color_GN84, color_GO4A]


# %%
# import matplotlib.pyplot as plt
# %matplotlib inline

zone_nord = ['eurat', 'nordat', "hn20", "glob"]
expe = ['GL55', 'GN84', 'GO4A']
color = [color_GL55, color_GN84, color_GO4A]

fig, axs = plt.subplots(nrows=1,  ncols=4,
                        sharex='col', sharey='row',
                        figsize=(35, 5), dpi=300)  # 46,5


for izone, zon in enumerate(zone):
    print(izone, zon)
    for isim, sim in enumerate(expe):
        print(isim, sim)

        df_ech = pd.read_csv(Dir+str(sim)+'_'+str(zon) +
                             '_prevision_bouees_flag_ech.csv')
        df_ech.head()
        df_ech.dtypes

        df_ech = df_ech.loc[((df_ech['expe'] == sim))]

        df_ech.drop(index=df_ech.index[0], axis=0, inplace=True)
        # fct get only one time repeat variables
        ech2 = np.unique(df_ech['ech'])
        N_level = ['SFX.TEMP_OC1']

        # Make mean over each forecast term - Biais
        data = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
                'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
                'biais_obs': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_res = pd.DataFrame(data)  # Create an empty DF

        i = 0
        for e in ech2:
            for N in N_level:
                biais_mean = float(df_ech.loc[(df_ech['ech'] == e) & (
                    df_ech['level'] == -0.494025)].mean()['biais_obs'])
                # print(biais_mean)
                df_res.iloc[i] = np.array([N, e, round(biais_mean, 5)])
                i = i + 1

        # Make mean over each forecast term - Buyos Nb
        data2 = {'level': np.repeat(np.nan, len(ech2)*len(N_level)),
                 'ech': np.repeat(np.nan, len(ech2)*len(N_level)),
                 'Nb_obs': np.repeat(np.nan, len(ech2)*len(N_level))}
        df_buoys = pd.DataFrame(data2)

        i = 0
        for e in ech2:
            for N in N_level:
                Nb_bouees_mean = float(df_ech.loc[(df_ech['ech'] == e) & (
                    df_ech['level'] == -0.494025)].mean()['Nb_bouees_an'])
                # print(biais_mean)
                df_buoys.iloc[i] = np.array([N, e, round(Nb_bouees_mean, 3)])

        print('min of biais', df_res.biais_obs.min())
        min = float(df_res.biais_obs.min())
        print('max of biais', df_res.biais_obs.max())
        max = float(df_res.biais_obs.max())

        # print('min of buoys numbers',df_buoys.Nb_obs.min())
        # min=float(df_buoys.Nb_obs.min())
        # print('max of buoys numbers',df_buoys.Nb_obs.max())
        # max=float(df_buoys.Nb_obs.max())

        df_res["biais_obs"] = df_res["biais_obs"].astype(float)
        df_buoys["Nb_obs"] = df_buoys["Nb_obs"].astype(float)

        # while izone <= 1:
        #     ax = axs[izone*0][izone]  # names axes for each zone
        # else:
        #     ax = axs[1][izone-2]
        
        ax=axs[izone]

        ax.plot(df_res.ech, df_res.biais_obs, color=color[isim], linewidth=2)
        ax.set_xticks(df_res.ech[np.arange(17)],
                      labels=np.arange(0, 102, 6))

        start = -0.05
        end = 0.7
        stepsize = 0.1
        ax.yaxis.set_ticks(np.arange(start, end, stepsize))

        axs[0].set_ylabel('Biais °C', fontsize=26)
        ax.set_xlabel('Time (h)', fontsize=26)
        ax.tick_params(axis='y', labelsize=14)
        ax.tick_params(axis='x', labelsize=14)

        ax.set_title('Zone:'+str(zon), fontsize=28, y=1.05)
        
        
         
        #diff_Number of buoys
        ax2 = axs[izone].twinx()
        ax2.plot(df_res.ech,df_buoys.Nb_obs,color=color[isim], linestyle=':', alpha=1)
        stepsize2=100
        start2=50
        end2=600
        ax2.yaxis.set_ticks(np.arange(start2, end2, stepsize2))
        ax2.tick_params(axis='y', labelsize=14)

        fig.suptitle(
            'Two months mean biais OML mod. against in situ data', fontsize=28, y=1.15)
        fig.legend(['L.26', 'Xp.finale', 'Xp.buyos'],
                   fontsize=30, bbox_to_anchor=(1.1, 1.05))
        plt.legend()
        fig.savefig(Dir_im+'biais_bouees_allzones.png', bbox_inches="tight")
