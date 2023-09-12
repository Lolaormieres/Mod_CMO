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

cterm='0'
cterm_for='96'
## Definition color
color_GL55='#648FFF'  ## GL55
color_GKOR='#785EF0'  ## GKOR
color_GKPH='#DC267F'  ## GKPH
color_GM6E='#FE6100'  ## GM6E
color_GMOT='grey'  ## GM6E

# =============================================================================
# ### Read DF an-mer
# =============================================================================
#csv_file = tf.keras.utils.get_file('GMOT/septembre_GMOT_an-mer_Nord_Atlantique.csv')
zone = 'nord-atlantique'
nom_zone = 'Nord-Atlantique'
df_an = pd.read_csv('DF_3/septembre_octobre__an-mer_'+str(zone)+'.csv')
df_an.head()
df_an.dtypes
# =============================================================================
# ## Analyse - Mercator
# =============================================================================
niv='SFX.TEMP_OC12'

## Vmin, Vmax
print('min tot',df_an.biais.min())
print('max tot',df_an.biais.max())
vmin=df_an.biais.min()
vmax=df_an.biais.max()

if abs(vmin)>abs(vmax):
    vmax=abs(vmin)
else:
    vmin=-vmax

# xp=['GL55','GKOR','GM6E','GMOT','GKPH']
xp=['GKPH','GKOR']
sim=['GKPH','GKOR']#],'GKPH']
nom_sim=['R.50','L.12']

# for isim, sim in enumerate(xp):
# print(isim,sim)
fig, ax = plt.subplots(2, 1, figsize=(8, 7),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig.tight_layout(pad=3)  

for isim, sim in enumerate(xp):
    print(isim,sim)
    ### TEST
    %matplotlib inline
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-white')
    import numpy as np
    N=5
    
    ### DATE - pour un niveau => 29 jours du 02 au 30
    df_an = pd.read_csv('DF_3/septembre_octobre__an-mer_'+str(zone)+'.csv')
    df_an.head()
    df_an.dtypes
    
    df_an_sim=df_an.loc[((df_an['expe']==sim) )][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
    #print(df_GL55)
    
    
    ## Prends Dates GKOR pour x labels
    df_date_labels=df_an.loc[((df_an['expe']=='GKOR') )][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
    #print(df_GL55)
    df_date_labels=df_date_labels.loc[((df_date_labels['var']=='SFX.TEMP_OC1') )][['date']]
    
    # Convert the date to datetime64
    df_an_sim['date'] = pd.to_datetime(df_an_sim['date'], format='%d/%m/%y')
    mask = (df_an_sim['date'] < '19/10/2022')  
    print(df_an_sim.loc[mask])
    
    ## DATE GKOR mask
    
    df_an_sim_corr=df_an_sim.loc[mask]
    
    import numpy as np
    import pandas as pd
    df_an_date=df_an_sim_corr.loc[((df_an_sim_corr['var']=='SFX.TEMP_OC1') )][['date']]
    
    
    # df_an_GKOR=df_an.loc[((df_an['expe']=='GKOR') )][['date','expe','var','level','biais','eqm','ech','zone','compar']]
    # df_an_date=df_an_GKOR.loc[((df_an_GKOR['var']=='SFX.TEMP_OC1') )][['date']]
    
    # df_an_final=df_an_GKPH.loc[(df_an_GKPH['date']==df_an_date) ][['date','expe','var','level','biais','eqm','ech','zone','compar']]
    #print(df_an_date) 
    
    
    ## DATE - XTICK
    jours=len(df_an_date)
    ## xLabels - Dates
    L=len(df_an_date)
    lst = list(np.arange(1,L+1))
    xx=lst
    print('lst',lst)
    xind = xx[0:len(xx):N] # Pas des dates plot
    # print(xind)
    
    
    ## Liste levels (12 niveaux)
    df_an_level=df_an_sim_corr.loc[((df_an_sim_corr['date']=='2022-09-02') )][['level']]
    levels=[]
    levels=df_an_level
    levels=np.array(levels)
    
    ## PLOT CONTOUR
    x = lst
    y = levels
    X, Y = np.meshgrid(x, y)
    
    biais=[]
    biais=df_an_sim_corr.biais
    print('min',min(biais))
    print('max',max(biais))
    biais=np.array(biais)
    # print('arr biais',biais)
    # print(len(biais))
    
    ## Reshape bias en 2D
    reshaped_biais = biais.reshape(jours,12).T
    Z = reshaped_biais
    # print(Z)
    
    # Basic contour plot
    fig, ax = plt.subplots()
    xlabels = df_date_labels.date[0:len(df_date_labels.date):N]
    plt.xticks(xind,labels=xlabels,rotation=45,fontsize=14)
    plt.title("Biais de température entre l'analyse (exp: "+str(nom_sim[isim])+") \n et Mercator, "+str(nom_zone), fontsize=14)
    plt.xlabel('Date', fontsize=15)
    ax.yaxis.set_tick_params(labelsize=13)
    #ax[0].set.xlabel('Dates')
    plt.ylabel('Niveaux (m)', fontsize=15)
    #CS = ax.contourf(X, Y, Z,levels=[-2,-0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1,2],cmap=plt.cm.bwr)
    #CS1 = ax.contour(X, Y, Z,cmap=plt.cm.bwr, levels=[-0.3,-0.25,-0.2,-0.15,-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30])
    #CS1 = ax.contour(X, Y, Z,cmap='seismic', levels=[-0.45,-0.40,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1, -0.05,-0.01,0.01, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45])
    # CS1 = ax.contour(X, Y, Z,cmap='seismic', levels=np.linspace(round(vmin,2),round(vmax,2),19)) ## N ATLANTIQUE
    CS1 = ax.contourf(X, Y, Z,cmap='seismic', levels=[-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4]) ## N ATLANTIQUE
    # CS1 = ax.contour(X, Y, Z,cmap='seismic', levels=np.linspace(round(vmin,2),round(vmax,2),17)) ## TROPIQUES
    # CS1 = ax.contour(X, Y, Z,cmap='seismic', levels=np.linspace(round(vmin,1),round(vmax,1),17)) ## EURAT
    # CS_lab = ax.contour(X, Y, Z,cmap='seismic', levels=np.linspace(round(vmin,1),round(vmax,1),8))
    # manual_locations = [-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4]
    cbar=plt.colorbar(CS1)
    # cbar.ax.set_yticklabels(fontsize=15)
    # cbar.set_label('Biais (°C)', rotation=270)
    # plt.clabel(CS, inline=1, fontsize=10, manual=manual_locations
    cbar.set_label('Biais (°C)', fontsize=16)
    # cbar.ax.tick_params(labelsize=18)
    # cbar.ax.tick_params(labelsize='large')
    # cb=plt.colorbar(im,orientation='horizontal').set_label(label='Label',size=20,weight='bold')
    rmin=round(vmin)
    rmax=round(vmax)
    plt.clabel(CS1,inline=2, fontsize=16)
    plt.show()
    
    figname = '/cnrm/recyf/Data/users/ormieresl/SST_Diff_temporelle_an-mercator_'+str(nom_sim[isim])+'_'+str(cterm)+'UTC_'+str(zone)+'.png'
    fig.savefig(figname,dpi=300, format='png',bbox_inches='tight')
    figname = '/cnrm/recyf/Data/users/ormieresl/SST_Diff_temporelle_an-mercator_'+str(nom_sim[isim])+'_'+str(cterm)+'UTC_'+str(zone)+'.pdf'
    fig.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')
    
    # # Figure 1
    # # plt.figure()
    # xlabels = df_date_labels.date[0:len(df_date_labels.date):N]
    # ax[0].set_xticks(xind,labels=xlabels,rotation=45,fontsize=14)
    # #plt.xticks(rotation=45,fontsize=8)
    # # ax[0].set_yticks(fontsize=8)
    # ax[0].set_xlabel('Date')
    # ax[0].set_ylabel('Biais en °C')
    # #ax[0].set_title(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
    # plt.suptitle("Biais de température entre l'analyse("+str(cterm)+'h) (exp: '+str(nom_sim[isim])+') et Mercator, '+str(nom_zone), fontsize=14)
    # ax[0].plot(df_GL55.date,df_GL55.biais,label='GL55',color=color_GL55)
    # ax[0].plot(df_GKOR.date,df_GKOR.biais,label='GKOR',color=color_GKOR)
    # ax[0].plot(df_GKPH.date,df_GKPH.biais,label='GKPH',color=color_GKPH)
    # ax[0].plot(df_GM6E.date,df_GM6E.biais,label='GM6E',color=color_GM6E)
    # ax[0].plot(df_GMOT.date,df_GMOT.biais,label='GMOT',color=color_GMOT)
    # ax[0].axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    # ax[0].legend(['GL55','GKOR','GKPH','GM6E','GMOT'])
    # # plt.show()
    # # plt.legend()
    # # figname = '/'+str(niv)+'_Diff_temporelle_prev-an_'+'_'+str(cterm)+'UTC'+str(zone)
    # # plt.savefig(figname, dpi=100, format='png',bbox_inches='tight')
    
    
    # # Figure 2
    # # plt.figure()
    # # Figure 1
    # # plt.figure()
    # xlabels = df_date_labels.date[0:len(df_date_labels.date):N]
    # ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=14)
    # #plt.xticks(rotation=45,fontsize=8)
    # # ax[0].set_yticks(fontsize=8)
    # ax[1].set_xlabel('Date')
    # ax[1].set_ylabel('Std en °C')
    # #ax[1].set_title(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
    # #plt.suptitle(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
    # plt.suptitle(str(niv)+', Prevision_'+str(cterm_for)+'h(xp) - Mercator, '+str(zone))
    # ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55)
    # ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR)
    # ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH)
    # ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E)
    # ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT)
    # ax[1].legend(['GL55','GKOR','GKPH','GM6E','GMOT'])
    # plt.show()



#%%
# =============================================================================
# READ DF PREV - MER
# =============================================================================
zone = 'nord-atlantique'
nom_zone = 'Nord-Atlantique'
### DATE - pour un niveau => 29 jours du 02 au 30
df_prev = pd.read_csv('DF_3/septembre_octobre__prev-mer_'+str(zone)+'.csv')
df_prev.head()
df_prev.dtypes

# =============================================================================
# ## Prevision - Mercator
# =============================================================================
niv='SFX.TEMP_OC12'

## Vmin, Vmax
print('min tot',df_an.biais.min())
print('max tot',df_an.biais.max())
vmin=df_prev.biais.min()
vmax=df_prev.biais.max()

if abs(vmin)>abs(vmax):
    vmax=abs(vmin)
else:
    vmin=-vmax

# xp=['GL55','GKOR','GM6E','GMOT','GKPH']
sim=['GKPH']#],'GKPH']
nom_sim='R.50'

xp=['GKPH','GKOR']
sim=['GKPH','GKOR']#],'GKPH']
nom_sim=['R.50','L.12']
# for isim, sim in enumerate(xp):
# print(isim,sim)


# print(isim,sim)
fig3, ax = plt.subplots(2, 1, figsize=(8, 7),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig3.tight_layout(pad=3)  
for isim, sim in enumerate(xp):
    # for isim, sim in enumerate(xp):
    
    print(isim,sim)
    ### TEST
    %matplotlib inline
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-white')
    import numpy as np
    N=5
    
    ### DATE - pour un niveau => 29 jours du 02 au 30
    df_prev = pd.read_csv('DF_3/septembre_octobre__prev-mer_'+str(zone)+'.csv')
    df_prev.head()
    df_prev.dtypes
    
    df_prev_sim=df_prev.loc[((df_prev['expe']==sim) )][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
    #print(df_GL55)
    
    
    ## Prends Dates GKOR pour x labels
    df_date_labels=df_prev.loc[((df_prev['expe']=='GKOR') )][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
    #print(df_GL55)
    df_date_labels=df_date_labels.loc[((df_date_labels['var']=='SFX.TEMP_OC1') )][['date']]
    
    # Convert the date to datetime64
    df_prev_sim['date'] = pd.to_datetime(df_prev_sim['date'], format='%d/%m/%y')
    mask = (df_prev_sim['date'] < '15/10/2022')  
    print(df_prev_sim.loc[mask])
    
    ## DATE GKOR mask
    
    df_prev_sim_corr=df_prev_sim.loc[mask]
    
    import numpy as np
    import pandas as pd
    df_prev_date=df_prev_sim_corr.loc[((df_prev_sim_corr['var']=='SFX.TEMP_OC1') )][['date']]
    
    
    # df_an_GKOR=df_an.loc[((df_an['expe']=='GKOR') )][['date','expe','var','level','biais','eqm','ech','zone','compar']]
    # df_an_date=df_an_GKOR.loc[((df_an_GKOR['var']=='SFX.TEMP_OC1') )][['date']]
    
    # df_an_final=df_an_GKPH.loc[(df_an_GKPH['date']==df_an_date) ][['date','expe','var','level','biais','eqm','ech','zone','compar']]
    #print(df_an_date) 
    
    
    ## DATE - XTICK
    jours=len(df_prev_date)
    ## xLabels - Dates
    L=len(df_prev_date)
    lst = list(np.arange(1,L+1))
    xx=lst
    print('lst',lst)
    xind = xx[0:len(xx):N] # Pas des dates plot
    # print(xind)
    
    
    ## Liste levels (12 niveaux)
    df_prev_level=df_prev_sim_corr.loc[((df_prev_sim_corr['date']=='2022-09-02') )][['level']]
    levels=[]
    levels=df_prev_level
    levels=np.array(levels)
    
    ## PLOT CONTOUR
    x = lst
    y = levels
    X, Y = np.meshgrid(x, y)
    
    biais=[]
    biais=df_prev_sim_corr.biais
    print('min',min(biais))
    print('max',max(biais))
    biais=np.array(biais)
    # print('arr biais',biais)
    # print(len(biais))
    
    ## Reshape bias en 2D
    reshaped_biais = biais.reshape(jours,12).T
    Z = reshaped_biais
    # print(Z)
    
    
    # Basic contour plot
    fig3, ax = plt.subplots()
    xlabels = df_date_labels.date[0:len(df_date_labels.date):N]
    plt.xticks(xind,labels=xlabels,rotation=45,fontsize=14)
    
    # plt.title("This is title number: " + r"$\bf{" + str(number) + "}$")
    
    plt.title("Biais de température prévue ("+str(cterm_for)+'h) (exp: '+str(nom_sim[isim])+', Rappel à 10%) \n par rapport à la REF.Mercator, Bassin '+str(nom_zone),fontsize=14)
    # plt.title("Biais de température prévue ("+str(cterm_for)+'h) (exp: '+str(nom_sim[isim])+', Rappel à 10%) \n par rapport à la REF.Mercator, Bassin'+str(nom_zone),fontsize=14,family="arial", backgroundcolor='lightgrey')
    plt.xlabel('Date',fontsize=14)
    #ax[0].set.xlabel('Dates')
    plt.ylabel('Niveaux (m)',fontsize=14)
    #CS = ax.contourf(X, Y, Z,levels=[-2,-0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1,2],cmap=plt.cm.bwr)
    #CS1 = ax.contour(X, Y, Z,cmap=plt.cm.bwr, levels=[-0.3,-0.25,-0.2,-0.15,-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30])
    #CS2 = ax.contour(X, Y, Z,cmap='bwr', levels=[-0.3,-0.25,-0.2,-0.15,-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30])
    echelon=np.arange(-0.7,0.7,0.5)
    #CS2 = ax.contour(X, Y, Z,cmap='seismic', levels=[-0.6,-0.55,-0.5,-0.45,-0.4, -0.35,-0.3,-0.25,-0.2,-0.15,-0.1, -0.05 , 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6])
    #CS2 = ax.contour(X, Y, Z,cmap='seismic', levels=np.linspace(round(vmin,1),round(vmax,1),15)) ## N ATLANTIQUE
    # levels=[-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4,0.5,0.6])
    CS2 = ax.contourf(X, Y, Z,cmap='seismic', levels=[-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6]) ## N ATLANTIQUE
    # CS2 = ax.contour(X, Y, Z,cmap='seismic', levels=np.linspace(round(vmin,1),round(vmax,1),9)) ## TROPIQUES
    # CS2 = ax.contour(X, Y, Z,cmap='seismic', levels=np.linspace(round(vmin,1),round(vmax,1),15))
    cbar=plt.colorbar(CS2)
    # cbar.set_label('Biais (°C)', rotation=270)
    cbar.set_label('Biais (°C)')
    plt.clabel(CS2,inline=2, fontsize=14)
    plt.show()
 
    
    
    figname = '/cnrm/recyf/Data/users/ormieresl/SST_Diff_temporelle_prevision-mercator_'+str(nom_sim[isim])+'_'+str(cterm_for)+'UTC_'+str(zone)+'.png'
    fig3.savefig(figname,dpi=300, format='png',bbox_inches='tight')
    figname = '/cnrm/recyf/Data/users/ormieresl/SST_Diff_temporelle_prevision-mercator_'+str(nom_sim[isim])+'_'+str(cterm_for)+'UTC_'+str(zone)+'.pdf'
    fig3.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')

    

    # ## PLOT propre
    
    ## TEST SUBPLOTS
    # print(isim,sim)
    # fig, ax = plt.subplots(2, 1, figsize=(8, 7),dpi=300)
    # # fig, ax = plt.subplots(2, 1)
    # fig.tight_layout(pad=3)  
    
    # ax[0].set_xticks(xind,labels=xlabels,rotation=45,fontsize=8)
    # #plt.xticks(rotation=45,fontsize=8)
    # # ax[0].set_yticks(fontsize=8)
    # ax[0].set_xlabel('Date')
    # ax[0].set_ylabel('Niveaux (m)',fontsize=14)
    # #ax[0].set_title(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
    # plt.suptitle("Biais de température entre la prévision ("+str(cterm_for)+'h) (exp: '+str(nom_sim[isim])+') et Mercator, '+str(nom_zone),fontsize=14)
    # CS2=ax[0].contourf(X, Y, Z,cmap='seismic', levels=[-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6]) 

    # cbar=plt.colorbar(CS2)
    # # # cbar.set_label('Biais (°C)', rotation=270)
    # # cbar.set_label('Biais (°C)')
    # plt.contourf(CS2,cmap='seismic')
    # plt.clabel(CS2,inline=2, fontsize=10)
    # # plt.show()
    # plt.legend()
    # # # figname = '/'+str(niv)+'_Diff_temporelle_prev-an_'+'_'+str(cterm)+'UTC'+str(zone)
    # # plt.savefig(figname, dpi=100, format='png',bbox_inches='tight')
    
    
    # # Figure 2
    # ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=8)
    # #plt.xticks(rotation=45,fontsize=8)
    # # ax[0].set_yticks(fontsize=8)
    # ax[1].set_xlabel('Date')
    # ax[1].set_ylabel('Std en °C')
    # #ax[1].set_title(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
    # #plt.suptitle(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
    # plt.suptitle(str(niv)+', Prevision_'+str(cterm_for)+'h(xp) - Mercator, '+str(zone))
    # ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55)
    # ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR)
    # ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH)
    # ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E)
    # ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT)
    # ax[1].legend(['GL55','GKOR','GKPH','GM6E','GMOT'])
    # plt.show()
    
    



























