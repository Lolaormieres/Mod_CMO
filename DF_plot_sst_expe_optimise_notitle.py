#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 14:34:45 2023

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
import matplotlib.axes as ax
import cartopy.crs as ccrs
import cartopy
import epygram
from datetime import datetime, timedelta
import vortex
import pandas as pd
import os
import numpy.ma as ma

import pandas as pd
epygram.init_env()


import matplotlib
for name, hex in matplotlib.colors.cnames.items():
    print(name, hex)


# ## Definition color
# color_GL55='#648FFF'  ## GL55
# color_GKOR='#785EF0'  ## GKOR
# color_GKPH='crimson'  ## GKPH
# color_GM6E='#FE6100'  ## GM6E
# color_GMOT='limegreen'  ## GMOT
# color_OPER='r'
# color_GK7C='black' 
# color_GN3C='gold'


## Definition color
color_GL55='royalblue' ## GL55 '#648FFF'
color_GKOR='#785EF0'  ## GKOR
color_GKPH= 'chocolate' ## GKPH '#DC267F'
color_GM6E='#FE6100'  ## GM6E
color_GMOT='limegreen'  ## GMOT '#009E73'
color_GK7C='black' 
color_GKCJ='dimgrey' 
color_OPER='r'
color_GN3C='gold'
color_GN84='red'
color_GNBA='royalblue'
color_GO4A='orange'
color_GNYG='limegreen'
color_GNSR='limegreen'


## Definition color Pale
color_GL55='steelblue' ## GL55 '#648FFF'
color_GKOR='hotpink'  ## GKOR
color_GKPH= 'mediumpurple' ## GKPH '#DC267F'
color_GM6E='paleturquoise'
color_GMOT='peachpuff'#'bisque'#'peachpuff'
color_GN84= 'palegreen'#'red'
color_GO4A='rosybrown'
color_GNSR='palegreen'
color_GNYG='paleyellow'
color_GOJQ='navajowhite'


cterm='0'
cterm_for='96'
# niv='SFX.SST'
niv='SFX.SST'

dir_df='/home/ormieresl/Routines/DF_SST_mapfactor_mask/'
dir_fig='/cnrm/recyf/Data/users/ormieresl/plot_evol_sst_technical/'
dir_fig='/cnrm/recyf/Data/users/ormieresl/plot_evol_sst_soutenance/plot_evol_sst_soutenance_2/'
#%%

# =============================================================================
# ### Read DF an-mer
# =============================================================================
#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [0,360,-20,20],"hn20":[0,360,20,90],"hs20":[0,360,-90,-20],"glob":[0,360,-90,90],"med": [-3,16,30,44]}
compar='an-mer'

zone=[]


### Ceate Dataframe
raw_data = {'zone': [''],
                'expe': [''],
                'biais': [''],
                'eqm': ['']}
df = pd.DataFrame(raw_data, columns = ['zone','expe','biais','eqm'])


for key in zones:
  
    print(key)
    zone.append(key)
print(zone)
    
# LOOP XP, mean due date
for izone, zone in enumerate(zone):
    print(izone, zone)
    
    df_GL55 = pd.read_csv(dir_df+'GL55_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GL55.head()
    df_GL55.dtypes
    df_GL55.drop(index=df_GL55.index[0], axis=0, inplace=True)
    
    df_GKPH = pd.read_csv(dir_df+'GKPH_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GKPH.head()
    df_GKPH.dtypes
    df_GKPH.drop(index=df_GKPH.index[0], axis=0, inplace=True)
    
    df_GKOR = pd.read_csv(dir_df+'GKOR_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GKOR.head()
    df_GKOR.dtypes
    df_GKOR.drop(index=df_GKOR.index[0], axis=0, inplace=True)
    
    df_GM6E = pd.read_csv(dir_df+'GM6E_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GM6E.head()
    df_GM6E.dtypes
    df_GM6E.drop(index=df_GM6E.index[0], axis=0, inplace=True)
    
        
    df_GMOT = pd.read_csv(dir_df+'GMOT_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GMOT.head()
    df_GMOT.dtypes
    df_GMOT.drop(index=df_GMOT.index[0], axis=0, inplace=True)
    
    df_GN84 = pd.read_csv(dir_df+'GN84_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GN84.head()
    df_GN84.dtypes
    df_GN84.drop(index=df_GN84.index[0], axis=0, inplace=True)
    
    df_GO4A = pd.read_csv(dir_df+'GO4A_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GO4A.head()
    df_GO4A.dtypes
    df_GO4A.drop(index=df_GO4A.index[0], axis=0, inplace=True)
    
    df_GNSR = pd.read_csv(dir_df+'GNSR_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GNSR.head()
    df_GNSR.dtypes
    df_GNSR.drop(index=df_GNSR.index[0], axis=0, inplace=True)
    
    # #ONLY prev-mer??
    # df_GNYG = pd.read_csv(dir_df+'GNYG_'+str(zone)+'_'+str(compar)+'_sst.csv')
    # df_GNYG.head()
    # df_GNYG.dtypes
    # df_GNYG.drop(index=df_GNYG.index[0], axis=0, inplace=True)
    

    df_GKCJ = pd.read_csv(dir_df+'GKCJ_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GKCJ.head()
    df_GKCJ.dtypes
    df_GKCJ.drop(index=df_GKCJ.index[0], axis=0, inplace=True)
    
        
    df_GOJQ = pd.read_csv(dir_df+'GOJQ_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GOJQ.head()
    df_GOJQ.dtypes
    df_GOJQ.drop(index=df_GOJQ.index[0], axis=0, inplace=True)
    
    # FIGURE Save CNRM
    # %matplotlib inline
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-white')
    import numpy as np
    # N=2
    
    fig1, ax = plt.subplots(2, 1, figsize=(14, 9),dpi=300)
    # fig, ax = plt.subplots(2, 1)
    fig1.tight_layout(pad=28)  
    
    # Plot dates axe x
    N=12
    
    L=len(df_GNSR.date)
    print(L)
    lst = list(np.arange(1,L+1))
    print(lst)
    xx=lst
    xind = xx[0:len(xx):N] # Pas des dates plot
    xlabels = df_GNSR.date[0:len(df_GNSR.date):N]
    
    ## ARRANGE FIGURE
    ## BIAIS
    ax[0].set_xticks([])
    ax[0].set_ylabel('Bias (°C)',fontsize=24)
    ax[0].yaxis.set_tick_params(labelsize=24)
    
    # plt.suptitle('SST ARPEGE analysée - SST Mercator (REF) \n\n',fontsize=18)
    # plt.suptitle('Analysed SST - Instantaneous Mercator SST (REF) \n zone: '+str(nom_zone),fontsize=26, y=0.96)
    # plt.suptitle('Analysed SST - Instantaneous Mercator SST (REF) \n'+str(zone),fontsize=26, y=0.97)
    
    ax[0].plot(df_GL55.date,df_GL55.biais,label='L.26',color=color_GL55,linewidth=3)
    ax[0].plot(df_GKOR.date,df_GKOR.biais,label='L.12',color=color_GKOR,linewidth=3)
    ax[0].plot(df_GKPH.date,df_GKPH.biais,label='R.50',color=color_GKPH,linewidth=3)
    # ax[0].plot(df_GM6E.date,df_GM6E.biais,label='no.curr',color=color_GM6E,linewidth=3)
    # ax[0].plot(df_GMOT.date,df_GMOT.biais,label='min.bathy',color=color_GMOT,linewidth=3)
    # ax[0].plot(df_GN84.date,df_GN84.biais,label='Xp.OML',color=color_GN84,linewidth=3)
    # ax[0].plot(df_GO4A.date,df_GO4A.biais,label='GO4A',color=color_GO4A,linewidth=3)
    # ax[0].plot(df_GNSR.date,df_GNSR.biais,label='GNSR',color=color_GNSR,linewidth=3)
    # ax[0].plot(df_GOJQ.date,df_GOJQ.biais,label='GOJQ',color=color_GOJQ,linewidth=3)
    # ax[0].plot(df_GNYG.date,df_GNYG.biais,label='GNYG',color=color_GNYG,linewidth=3)
    # ax[0].plot(df_OPER.date,df_OPER.biais,label='OPER.OSTIA',color=color_OPER)
    ax[0].plot(df_GKCJ.date,df_GKCJ.biais,label='Double GKCJ',color=color_GKCJ,linewidth=3)


    ax[0].axhline(y=0, color='black', linestyle='--', linewidth=2)
    # ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)
        
    #Nb points - stretch grid
    ax2 = ax[0].twinx()
    ax2.set_ylim(0,3000000)
    ax2.plot(df_GN84.date,df_GN84.Nb_data_diff,color='lightgrey',linestyle=':',linewidth=3)
    ax2.set_ylabel('Grid points (10^6)',color='grey', fontsize=24)
    ax2.tick_params(axis="y", labelcolor='dimgrey',labelsize=24)
    
    
    ## STD
    # xlabels = df_GNSR.date[0:len(df_GNSR.date):N]
    ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=22)
    # ax[1].set_xlabel('Date',fontsize=22)
    ax[1].set_ylabel('Standard deviation (°C)',fontsize=24)
    ax[1].yaxis.set_tick_params(labelsize=24)

    ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55,linewidth=3)
    ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR,linewidth=3)
    ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH,linewidth=3)
    # ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E,linewidth=3)
    # ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT,linewidth=3)
    # ax[1].plot(df_GN84.date,df_GN84.eqm,label='Xp.OML',color=color_GN84,linewidth=3)
    # ax[1].plot(df_GO4A.date,df_GO4A.eqm,label='GO4A.buyos',color=color_GO4A,linewidth=3)
    # ax[1].plot(df_GNSR.date,df_GNSR.eqm,label='GNSR.wind',color=color_GNSR,linewidth=3)
    # ax[1].plot(df_GNYG.date,df_GNYG.eqm,label='GNYG.wind.R10',color=color_GNYG,linewidth=3)
    # ax[1].plot(df_GOJQ.date,df_GOJQ.eqm,label='GOJQ.R10x4',color=color_GOJQ,linewidth=3)
    # ax[1].plot(df_OPER.date,df_OPER.eqm,label='OPER.OSTIA',color=color_OPER)
    ax[1].plot(df_GKCJ.date,df_GKCJ.eqm,label='Double GKCJ',color=color_GKCJ,linewidth=3)

    # ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.OMLe','SST.cst.mer'],fontsize=14, loc='upper right')
    #ax[1].legend(['L.26','L.12','R.50','Xp.OMLe','GO4A','GNSR.wind','oper.2024'],fontsize=19, loc='best')
    #ax[1].legend(['L.26','L.12','R.50','L.12.noccurr','L.12.nocurr.drown','Xp.OMLe','Xp.buyos','Xp.wind','Xp.wind.R10','oper.2024'],fontsize=19, loc='best')
    # ax[1].legend(['L.26','L.12','R.50','L.12.noccurr','L.12.nocurr.drown','Xp.OMLe','Xp.buyos','Xp.wind','Xp.R10x4','oper.2024'],fontsize=19,bbox_to_anchor=(1.1, 1.05))
    # ax[1].legend(['L.26','L.12','R.50','Xp.OMLe','Xp.buyos','Xp.R10x4','Ref.noOML'],fontsize=19,bbox_to_anchor=(1.1, 1.05))
    ax[1].legend(['L.26','L.12','R.50','Xp.OML','Ref.noOML'],fontsize=27,bbox_to_anchor=(1.02, 1.01))

    ax[1].set_ylim(0, 0.7)
    
    biais_moy=np.round(df_GKCJ[97:159].biais.mean(),1)
    std_moy=np.round(df_GKCJ[97:159].eqm.mean(),1)
    new_row = {'zone':zone,'expe':'GKCJ','biais':biais_moy, 'eqm':std_moy}    
    df.loc[len(df)] = new_row
      

    # SUPRIMER 1er COL DF
    #df.drop(index=df_diurne_an.index[0], axis=0, inplace=True) ## DELETE FIR

    df.to_csv(dir_fig+'analyse_sst_12dec-12fev_REF.csv')
    
    
    plt.show()
    figname = dir_fig+'SST_analyse-mercator'+'_'+str(zone)+'.png'
    fig1.savefig(figname,dpi=150, format='png',bbox_inches='tight')
    figname = dir_fig+'SST_analyse-mercator'+'_'+str(zone)+'.pdf'
    fig1.savefig(figname,dpi=150, format='pdf',bbox_inches='tight')


#%%

# =============================================================================
# ### Read DF an-mer
# =============================================================================
#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [0,360,-20,20],"hn20":[0,360,20,90],"hs20":[0,360,-90,-20],"glob":[0,360,-90,90],"med": [-3,16,30,44]}
compar='prev-mer'
zone=[]
### Ceate Dataframe
raw_data = {'zone': [''],
                'expe': [''],
                'biais': [''],
                'eqm': [''],
                'Nb': ['']}
df2 = pd.DataFrame(raw_data, columns = ['zone','expe','biais','eqm','Nb'])


for key in zones:
    print(key)
    zone.append(key)
print(zone)

# LOOP XP, mean due date
for izone, zone in enumerate(zone):
    print(izone, zone)
    df_GL55 = pd.read_csv(dir_df+'GL55_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GL55.head()
    df_GL55.dtypes
    df_GL55.drop(index=df_GL55.index[0], axis=0, inplace=True)
    
    df_GKPH = pd.read_csv(dir_df+'GKPH_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GKPH.head()
    df_GKPH.dtypes
    df_GKPH.drop(index=df_GKPH.index[0], axis=0, inplace=True)
    
    df_GKOR = pd.read_csv(dir_df+'GKOR_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GKOR.head()
    df_GKOR.dtypes
    df_GKOR.drop(index=df_GKOR.index[0], axis=0, inplace=True)
    
    df_GN84 = pd.read_csv(dir_df+'GN84_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GN84.head()
    df_GN84.dtypes
    df_GN84.drop(index=df_GN84.index[0], axis=0, inplace=True)
    
    df_GM6E = pd.read_csv(dir_df+'GM6E_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GM6E.head()
    df_GM6E.dtypes
    df_GM6E.drop(index=df_GM6E.index[0], axis=0, inplace=True)
    
        
    df_GMOT = pd.read_csv(dir_df+'GMOT_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GMOT.head()
    df_GMOT.dtypes
    df_GMOT.drop(index=df_GMOT.index[0], axis=0, inplace=True)
    
    df_GO4A = pd.read_csv(dir_df+'GO4A_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GO4A.head()
    df_GO4A.dtypes
    df_GO4A.drop(index=df_GO4A.index[0], axis=0, inplace=True)
    
    
    df_GNSR = pd.read_csv(dir_df+'GNSR_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GNSR.head()
    df_GNSR.dtypes
    df_GNSR.drop(index=df_GNSR.index[0], axis=0, inplace=True)
    
    # df_GNYG = pd.read_csv(dir_df+'GNYG_'+str(zone)+'_'+str(compar)+'_sst.csv')
    # df_GNYG.head()
    # df_GNYG.dtypes
    # df_GNYG.drop(index=df_GNYG.index[0], axis=0, inplace=True)
     
    df_GKCJ = pd.read_csv(dir_df+'GKCJ_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GKCJ.head()
    df_GKCJ.dtypes
    df_GKCJ.drop(index=df_GKCJ.index[0], axis=0, inplace=True)
    
    df_GOJQ = pd.read_csv(dir_df+'GOJQ_'+str(zone)+'_'+str(compar)+'_sst.csv')
    df_GOJQ.head()
    df_GOJQ.dtypes
    df_GOJQ.drop(index=df_GOJQ.index[0], axis=0, inplace=True)
    
    # FIGURE Save CNRM
    %matplotlib inline
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-white')
    import numpy as np
    # N=2
    
    fig2, ax = plt.subplots(2, 1, figsize=(14, 9),dpi=300)
    # fig, ax = plt.subplots(2, 1)
    fig2.tight_layout(pad=28)  
    
    # Plot dates axe x
    N=12
    
    L=len(df_GN84.date)
    print(L)
    lst = list(np.arange(1,L+1))
    print(lst)
    xx=lst
    xind = xx[0:len(xx):N] # Pas des dates plot
    xlabels = df_GN84.date[0:len(df_GN84.date):N]
    
    ## ARRANGE FIGURE
    ## BIAIS
    ax[0].set_xticks([])
    ax[0].set_ylabel('Bias (°C)',fontsize=24)
    ax[0].yaxis.set_tick_params(labelsize=24)
    
    # plt.suptitle('96-h forecasted SST - Instantaneous Mercator SST (REF) \n'+str(zone),fontsize=22, y=0.95)

    ax[0].plot(df_GL55.date,df_GL55.biais,label='L.26',color=color_GL55,linewidth=3)
    ax[0].plot(df_GKOR.date,df_GKOR.biais,label='L.12',color=color_GKOR,linewidth=3)
    ax[0].plot(df_GKPH.date,df_GKPH.biais,label='R.50',color=color_GKPH,linewidth=3)
    # ax[0].plot(df_GM6E.date,df_GM6E.biais,label='no.curr',color=color_GM6E,linewidth=3)
    # ax[0].plot(df_GMOT.date,df_GMOT.biais,label='min.bathy',color=color_GMOT,linewidth=3)
    ax[0].plot(df_GN84.date,df_GN84.biais,label='Xp.OML',color=color_GN84,linewidth=3)
    # ax[0].plot(df_GO4A.date,df_GO4A.biais,label='GO4A',color=color_GO4A,linewidth=3)
    # ax[0].plot(df_GNSR.date,df_GNSR.biais,label='GNSR',color=color_GNSR,linewidth=3)
    # ax[0].plot(df_GNYG.date,df_GNYG.biais,label='GNYG',color=color_GNYG,linewidth=3)
    # ax[0].plot(df_GOJQ.date,df_GOJQ.biais,label='GOJQ',color=color_GOJQ,linewidth=3)
    # ax[0].plot(df_OPER.date,df_OPER.biais,label='OPER.OSTIA',color=color_OPER)
    ax[0].plot(df_GKCJ.date,df_GKCJ.biais,label='Double GKCJ',color=color_GKCJ,linewidth=3)
    ax[0].axhline(y=0, color='black', linestyle='--', linewidth=1.7)
    # ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)
    
    #Nb points - stretch grid
    ax2 = ax[0].twinx()
    ax2.set_ylim(0,3000000)

    ax2.plot(df_GN84.date,df_GN84.Nb_data_diff,color='lightgrey',linewidth=3,linestyle=':')
    ax2.set_ylabel('Grid points (10^6)',color='grey', fontsize=22)
    ax2.tick_params(axis="y", labelcolor='dimgrey',labelsize=22)
    
    # print(df_GN84.Nb_data_diff[0])
    ## STD
    xlabels = df_GN84.date[0:len(df_GN84.date):N]
    ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=22)
    # ax[1].set_xlabel('Date',fontsize=22)
    ax[1].set_ylabel('Standard deviation (°C)',fontsize=24)
    ax[1].yaxis.set_tick_params(labelsize=24)

    ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55,linewidth=3)
    ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR,linewidth=3)
    ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH,linewidth=3)
    # ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E,linewidth=3)
    # ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT,linewidth=3)
    ax[1].plot(df_GN84.date,df_GN84.eqm,label='Xp.OML',color=color_GN84,linewidth=3)
    # ax[1].plot(df_GO4A.date,df_GO4A.eqm,label='GO4A.buyos',color=color_GO4A,linewidth=3)
    # ax[1].plot(df_GNSR.date,df_GNSR.eqm,label='GNSR.wind',color=color_GNSR,linewidth=3)
    # ax[1].plot(df_GNYG.date,df_GNYG.eqm,label='GNYG.wind.R10',color=color_GNYG,linewidth=3)
    # ax[1].plot(df_GOJQ.date,df_GOJQ.eqm,label='GOJQ.R10x4',color=color_GOJQ,linewidth=3)
    # ax[1].plot(df_OPER.date,df_OPER.eqm,label='OPER.OSTIA',color=color_OPER)
    ax[1].plot(df_GKCJ.date,df_GKCJ.eqm,label='Double GKCJ',color=color_GKCJ,linewidth=3)
    # ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.OMLe','SST.cst.mer'],fontsize=14, loc='upper right')
    #ax[1].legend(['L.26','L.12','R.50','L.12.noccurr','L.12.nocurr.drown','Xp.OMLe','Xp.buyos','Xp.wind','Xp.wind.R10','oper.2024'],fontsize=19)
   
    # ax2.set_ylim(45000,2500000)   
    # ax[0].set_ylim(-0.3, 1.5)
    # ax[1].set_ylim(0, 1.4)
    # ax[1].legend(['L.26','L.12','R.50','L.12.noccurr','L.12.nocurr.drown','Xp.OMLe','Xp.buyos','Xp.wind','Xp.R10x4','oper.2024'],fontsize=19,bbox_to_anchor=(1.1, 1.05))
    # ax[1].legend(['L.26','L.12','R.50','Xp.OMLe','Xp.buyos','Xp.R10x4','REF.noOML'],fontsize=19,bbox_to_anchor=(1.1, 1.05))
    # ax[1].legend(['L.26','L.12','R.50','REF.noOML'],fontsize=27,bbox_to_anchor=(1.02, 0.97))
    
    biais_moy=np.round(df_GKCJ[97:159].biais.mean(),1)
    std_moy=np.round(df_GKCJ[97:159].eqm.mean(),1)
    Nb_moy=np.round(df_GKCJ[97:159].Nb_data_diff.mean())
    
    new_row = {'zone':zone,'expe':'GKCJ','biais':biais_moy, 'eqm':std_moy,'Nb':Nb_moy}    
    df2.loc[len(df2)] = new_row
      
    # SUPRIMER 1er COL DF
    #df2.drop(index=df2.index[0], axis=0, inplace=True) ## DELETE FIR
    df2.to_csv(dir_fig+'forecast_sst_12dec-12fev_REF.csv')
    
    
    plt.show()
    figname = dir_fig+'SST_forecast-mercator'+'_'+str(zone)+'_soutenance2.png'
    fig2.savefig(figname,dpi=150, format='png',bbox_inches='tight')
    figname = dir_fig+'SST_forecast-mercator'+'_'+str(zone)+'_soutenance2.pdf'
    fig2.savefig(figname,dpi=150, format='pdf',bbox_inches='tight')
