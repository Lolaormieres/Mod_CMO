#%%
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



## Definition color
# color_GL55='#648FFF'  ## GL55
# color_GKOR='#785EF0'  ## GKOR
# color_GKPH='#DC267F'  ## GKPH
# color_GM6E='#FE6100'  ## GM6E
# color_GMOT='#009E73'  ## GMOT



## Definition color
color_GL55='royalblue' ## GL55 '#648FFF'
color_GKOR='#785EF0'  ## GKOR
color_GKPH= '#FF0000' ## GKPH '#DC267F'
color_GM6E='#FE6100'  ## GM6E
color_GMOT='limegreen'  ## GMOT '#009E73'
color_GK7C='dimgrey' 
color_OPER='r'
color_GN3C='gold'
color_GN84='midnightblue'
color_GNBA='lightseagreen'


zone = 'nord-atlantique'
nom_zone='Nord-Atlantique'
zone = 'eurat'
nom_zone='eurat'
zone = 'tropiques'
nom_zone='Tropiques'
niv='SFX.TEMP_OC1'

#%%
# =============================================================================
# ### Read DF an-mer
# =============================================================================
#csv_file = tf.keras.utils.get_file('GMOT/septembre_GMOT_an-mer_Nord_Atlantique.csv')

compar='an-mer'
cterm='0'
cterm_for='96'


df_an = pd.read_csv('DF_optimise/septembre_octobre_'+str(compar)+'_'+str(zone)+'.csv')
df_an.head()
df_an.dtypes


# =============================================================================
# ## Analyse - Mercator
# =============================================================================
xp='GL55'
df_GL55=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
print(df_GL55)

xp='GKOR'
df_GKOR=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
print(df_GKOR)

xp='GKPH'
df_GKPH=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
print(df_GKPH)

xp='GM6E'
df_GM6E=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']]
print(df_GM6E)

xp='GMOT'
df_GMOT=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']]
print(df_GMOT)


xp='GN84'
df_GN84=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']]
print(df_GN84)


      
# Convert the date for GN84 because start 27/07/2022 --> 02/09/2022
df_GN84['date'] = pd.to_datetime(df_GN84['date'], format='%d/%m/%y') ## CONVERT DATE FORMAT
print(df_GN84)
mask = (df_GN84['date'] > '2022-09-01')
print(df_GN84.loc[mask])
df_GN84_corr=df_GN84.loc[mask] ## MASK SUR LES DATES

## RECONVERT TO %d/%m/%y
df_GN84_corr["date"] = pd.to_datetime(df_GN84_corr["date"],errors="coerce").dt.strftime("%d/%m/%y")
print(df_GN84_corr)
# ### TEST
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
# N=2

fig, ax = plt.subplots(2, 1, figsize=(8, 7),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig.tight_layout(pad=3)  

# Plot dates axe x
N=5

L=len(df_GL55.date)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
# Figure 1
# plt.figure()
xlabels = df_GL55.date[0:len(df_GL55.date):N]
ax[0].set_xticks(xind,labels=xlabels,rotation=45,fontsize=8)


#plt.xticks(rotation=45,fontsize=8)
# ax[0].set_yticks(fontsize=8)
ax[0].set_xlabel('Date')
ax[0].set_ylabel('Différence en °C')
# ax[0].set_title(str(niv)+'_A(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
# plt.suptitle(str(niv)+'_A(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
ax[0].plot(df_GL55.date,df_GL55.biais,label='GL55',color=color_GL55)
ax[0].plot(df_GKOR.date,df_GKOR.biais,label='GKOR',color=color_GKOR)
ax[0].plot(df_GKPH.date,df_GKPH.biais,label='GKPH',color=color_GKPH)
ax[0].plot(df_GM6E.date,df_GM6E.biais,label='GM6E',color=color_GM6E)
ax[0].plot(df_GMOT.date,df_GMOT.biais,label='GMOT',color=color_GMOT)
ax[0].plot(df_GN84_corr.date,df_GN84_corr.biais,label='GN84',color=color_GN84)
ax[0].axhline(y=0, color='black', linestyle='--', linewidth=0.5)
# ax[0].legend(['GL55','GKOR','GKPH','GM6E','GMOT'])
# plt.show()
# plt.legend()
# figname = '/'+str(niv)+'_Diff_temporelle_prev-an_'+'_'+str(cterm)+'UTC'+str(zone)
# plt.savefig(figname, dpi=100, format='png',bbox_inches='tight')


# Figure 2
# plt.figure()
# Figure 1
# plt.figure()
xlabels = df_GL55.date[0:len(df_GL55.date):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=8)
#plt.xticks(rotation=45,fontsize=8)
# ax[0].set_yticks(fontsize=8)
ax[1].set_xlabel('Date')
ax[1].set_ylabel('Std °C')
# ax[1].set_title(str(niv)+'_A(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
# plt.suptitle(str(niv)+'_A(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
plt.suptitle(str(niv)+', Analyse_'+str(cterm)+'h(xp) - Mercator, '+str(nom_zone))
ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55)
ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR)
ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH)
ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E)
ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT)
ax[1].plot(df_GN84_corr.date,df_GN84_corr.eqm,label='GN84',color=color_GN84)
ax[1].legend(['GL55','GKOR','GKPH','GM6E','GMOT','GN84'])
plt.show()
# plt.legend()
# figname = str(expe)+'/'+str(niv)+'_Posqter_Diff_temporelle_prev-an_'+'_'+str(cterm_for)+'UTC'+str(zone)
# plt.savefig(figname, dpi=100, format='png',bbox_inches='tight')
# print(df_niv.biais)


#%%
# =============================================================================
# ### Read DF an-mer
# =============================================================================
compar='prev-mer'
cterm='0'
cterm_for='96'
# zone = 'nord_atlantique'

df_an = pd.read_csv('DF_optimise/septembre_octobre_'+str(compar)+'_'+str(zone)+'.csv')
df_an.head()
df_an.dtypes


#%%
# =============================================================================
# ## Prevision - Mercator
# =============================================================================
xp='GL55'
df_GL55=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
print(df_GL55)

xp='GKOR'
df_GKOR=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
print(df_GKOR)

xp='GKPH'
df_GKPH=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']] #& (df_an['var']==niv)[['date','expe','var','level','biais','eqm','ech','zone','compar']]]
print(df_GKPH)

xp='GM6E'
df_GM6E=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']]
print(df_GM6E)

xp='GMOT'
df_GMOT=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']]
print(df_GMOT)


xp='GN84'
df_GN84=df_an.loc[((df_an['expe']==xp) & (df_an['var']==niv))][['date','expe','var','level','biais','eqm','ech','zone','compar']]
print(df_GN84)


      
# Convert the date for GN84 because start 27/07/2022 --> 02/09/2022
df_GN84['date'] = pd.to_datetime(df_GN84['date'], format='%d/%m/%y') ## CONVERT DATE FORMAT
print(df_GN84)
mask = (df_GN84['date'] > '2022-09-01')
print(df_GN84.loc[mask])
df_GN84_corr=df_GN84.loc[mask] ## MASK SUR LES DATES

## RECONVERT TO %d/%m/%y
df_GN84_corr["date"] = pd.to_datetime(df_GN84_corr["date"],errors="coerce").dt.strftime("%d/%m/%y")
print(df_GN84_corr)

    
# ### TEST
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
# N=2

fig, ax = plt.subplots(2, 1, figsize=(8, 7),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig.tight_layout(pad=3)  


# Plot dates axe x
N=5
L=len(df_GL55.date)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot

# Figure 1
# plt.figure()
xlabels = df_GL55.date[0:len(df_GL55.date):N]
ax[0].set_xticks(xind,labels=xlabels,rotation=45,fontsize=8)
#plt.xticks(rotation=45,fontsize=8)
# ax[0].set_yticks(fontsize=8)
ax[0].set_xlabel('Date')
ax[0].set_ylabel('Biais en °C')
#ax[0].set_title(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
plt.suptitle(str(niv)+' Prevision_'+str(cterm_for)+'(xp) - Mercator '+str(nom_zone))
ax[0].plot(df_GL55.date,df_GL55.biais,label='GL55',color=color_GL55)
ax[0].plot(df_GKOR.date,df_GKOR.biais,label='GKOR',color=color_GKOR)
ax[0].plot(df_GKPH.date,df_GKPH.biais,label='GKPH',color=color_GKPH)
ax[0].plot(df_GM6E.date,df_GM6E.biais,label='GM6E',color=color_GM6E)
ax[0].plot(df_GMOT.date,df_GMOT.biais,label='GMOT',color=color_GMOT)
ax[0].plot(df_GN84_corr.date,df_GN84_corr.biais,label='GN84',color=color_GN84)
ax[0].axhline(y=0, color='black', linestyle='--', linewidth=0.5)
# ax[0].legend(['GL55','GKOR','GKPH','GM6E','GMOT'])
# plt.show()
# plt.legend()
# figname = '/'+str(niv)+'_Diff_temporelle_prev-an_'+'_'+str(cterm)+'UTC'+str(zone)
# plt.savefig(figname, dpi=100, format='png',bbox_inches='tight')


# Figure 2
# plt.figure()
# Figure 1
# plt.figure()
xlabels = df_GL55.date[0:len(df_GL55.date):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=8)
#plt.xticks(rotation=45,fontsize=8)
# ax[0].set_yticks(fontsize=8)
ax[1].set_xlabel('Date')
ax[1].set_ylabel('Std en °C')
#ax[1].set_title(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
#plt.suptitle(str(niv)+'_P(xp) - A(Mer)_'+str(zone)+'_{}UTC'.format(cterm))
plt.suptitle(str(niv)+', Prevision_'+str(cterm_for)+'h(xp) - Mercator, '+str(nom_zone))
ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55)
ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR)
ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH)
ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E)
ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT)
ax[1].plot(df_GN84_corr.date,df_GN84_corr.eqm,label='GN84',color=color_GN84)
ax[1].legend(['GL55','GKOR','GKPH','GM6E','GMOT','GN84'])
plt.show()

# figname = str(expe)+'/'+str(niv)+'_Diff_temporelle_prev-an_'+str(expe)+'_'+str(cterm_for)+'UTC'+str(zone)
# plt.savefig(figname, dpi=100, format='png',bbox_inches='tight')
# print(df_niv.biais)

