#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 11:56:39 2023

@author: ormieresl
"""
#import bronx
from vortex import toolbox
import common, olive
import common.util.usepygram
import usevortex
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import epygram
from datetime import datetime,timedelta
import vortex
import os
import pandas as pd
import numpy.ma as ma
epygram.init_env()
#Pour traiter les fichiers au format netCDF
from netCDF4 import Dataset
#Imports pour tracer des cartes
import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.mpl.ticker as cticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Path to stock results
dir_df='/home/ormieresl/Routines/DF_bouees_ech_vf/'
dir_im='/cnrm/recyf/Data/users/ormieresl/plot_an_fg_depar_firstexp/'

expe=['GKCJ']
version='double'

# reference
expe_mer=['GO48']
compar='prev-mer'+expe_mer[0]

# Levels numbers
nb_level=1

# Start date
date_init = datetime(2022, 9,2, 0)
deltat = timedelta(days=1)

# Days numbers
nb=30
nbDate=nb # buyos dates
period = timedelta(days=nb)

# échéances
ech_fin=102
ech_ini=0
ech_step=6

# cterm_for='96'
zone =  'Glo'
period='sept'
# period='sept'
N=8
nbDate=120

# lecture des variables
name_sst=['SFX.SST']
name=['SFX.TEMP_OC'+str(i) for i in range(1,nb_level+1)]
name_sal=['SFX.SALT_OC'+str(i) for i in range(1,nb_level+1)]
name_mer=['SURF.THETAO'+str(i) for i in range(1,nb_level+1)]
name_sal_mer=['SURF.SALINO'+str(i) for i in range(1,nb_level+1)]
name_u=['SFX.UCUR_OC'+str(i) for i in range(1,nb_level+1)]
name_v=['SFX.VCUR_OC'+str(i) for i in range(1,nb_level+1)]
## Lecture des niveaux de la CMO
filename = 'level2.txt'
levelfic = np.loadtxt(filename, delimiter=',', dtype=float)
level = levelfic[0:nb_level]

## Definition color Pale
color_GL55='lightskyblue'#'paleturquoise' ## GL55 '#648FFF' 'steelblue'
color_GKOR='pink'#hotpink'  ## GKOR
color_GKPH= 'mediumpurple' ## GKPH '#DC267F'
color_GM6E='paleturquoise'
color_GMOT='peachpuff'
color_GN84='palegreen' # 'red'
color_GO4A='rosybrown'
color_GNSR='palegreen'
color_GNYG='yellow'
color_GOJQ='navajowhite'#burlywood' #'teal'
color_GKCJ='dimgrey'

#domaines lonmin,lonmax,latmin,matmax
zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [0,360,-20,20],"hn20":[0,360,20,90],"hs20":[0,360,-90,-20],"glob":[0,360,-90,90],"med": [-3,16,30,44]}
compar='an-mer'
#%%
## PLOT CYCLE DIURNE - BOUEES
expe=['GL55','GKOR','GKPH','GN84','GO4A','GOJQ','GKCJ']
XP='GN84'
SST_bouees=[]

### Ceate Dataframe
raw_data = {'date': [''],
        'expe': [''],
        'andepar': [''],
        'fgdepar': [''],
        'Nbbouees':[''],
        'Nbbouees_tot':[''],
        'Nbbouees_fg':[''],
        'Nbbouees_an':['']}
df_diurne_bouees= pd.DataFrame(raw_data, columns = ['date','expe','andepar','fgdepar','Nbbouees','Nbbouees_tot', 'Nbbouees_fg','Nbbouees_an'])

for e in range(len(expe)):
    print(e,expe[e])
    import datetime as dt
    year=2022
    month=9 
    day=2
    hour=0
    Delta=6
    Date_debut = dt.datetime(year, month, day, hour)
    timedelta = dt.timedelta(0,0,0,0,0,Delta)
    toto=[Date_debut+i*timedelta for i in range(nbDate)]
    # #Répertoire de travail
    chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(expe[e])+'/'
    
    #Loop: read files
    for date_dt in toto:
        print(date_dt)
        date = '%4.4i%2.2i%2.2iT%2.2i00A'%(date_dt.year,\
        date_dt.month,\
        date_dt.day,\
        date_dt.hour)
        print(date) #Read repertory name
        date_3 = '%2.2i%2.2i%2.2i%2.2i00'%(date_dt.year,\
        date_dt.month,\
        date_dt.day,\
        date_dt.hour)
        print(date_3)
        
        fic_nc = chemin_nc +'{}/surfan/odb-ecma.canari.surf.nc'.format(date)
        print(fic_nc)
        if not os.path.isfile(fic_nc):
            print('pas de fichier '+fic_nc)
    
        if os.path.isfile(fic_nc):
        #Ouverture du fichier netcdf au format dataset
            print('Ouverture du fichier',fic_nc)
            myfile = Dataset(fic_nc,'r')
        #print(myfile.variables) #structure fichier
    
        #Pour donner des noms explicites aux différentes variables :
        odb_key_val_dic = {}
        for variables in myfile.variables.keys():
        #instrument_int = myfile[variables].odb_name[:myfile[variables].odb_name.index('@')]
            odb_key_val_dic[variables] = myfile[variables][:]
                #print(variables,'->',instrument_int)
                
        # ind_flag = np.where((odb_key_val_dic['col_1']>1)&
        # (odb_key_val_dic['col_1']<2)&
        # (odb_key_val_dic['col_2']>40)&
        # (odb_key_val_dic['col_2']<42))
        
        #EURAT
        # zone='eurat'
        # ind_flag_tot = np.where((odb_key_val_dic['col_1']>-35)&
        # (odb_key_val_dic['col_1']<45)&
        # (odb_key_val_dic['col_2']>20)&
        # (odb_key_val_dic['col_2']<72))
        
        # ind_flag = np.where((odb_key_val_dic['col_1']>-35)&
        # (odb_key_val_dic['col_1']<45)&
        # (odb_key_val_dic['col_2']>20)&
        # (odb_key_val_dic['col_2']<72)&
        # (odb_key_val_dic['col_10']==1))
        
        # #Med
        # zone='med'
        # ind_flag_tot = np.where((odb_key_val_dic['col_1']>3)&
        # (odb_key_val_dic['col_1']<16)&
        # (odb_key_val_dic['col_2']>30)&
        # (odb_key_val_dic['col_2']<44))
        
        # ind_flag = np.where((odb_key_val_dic['col_1']>3)&
        # (odb_key_val_dic['col_1']<16)&
        # (odb_key_val_dic['col_2']>30)&
        # (odb_key_val_dic['col_2']<44)&
        # (odb_key_val_dic['col_10']==1))
        
        
        #N ATL
        # zone='nordat'
        # ind_flag_tot = np.where((odb_key_val_dic['col_1']<0)&
        # (odb_key_val_dic['col_1']>-80)&
        # (odb_key_val_dic['col_2']>20)&
        # (odb_key_val_dic['col_2']<70))
        
        # ind_flag = np.where((odb_key_val_dic['col_1']<0)&
        # (odb_key_val_dic['col_1']>-80)&
        # (odb_key_val_dic['col_2']>20)&
        # (odb_key_val_dic['col_2']<70)&
        # (odb_key_val_dic['col_10']==1))
        
        
        # TROPIQUES
        # zone='tropics'
        # ind_flag_tot = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']>-20)&
        # (odb_key_val_dic['col_2']<20))
        
        # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']>-20)&
        # (odb_key_val_dic['col_2']<20))
        # (odb_key_val_dic['col_10']==1))
        
                      
        ## South. Hem
        # zone='hs20'
        # ind_flag_tot = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']<-20)&
        # (odb_key_val_dic['col_2']>-75))
        
        # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']<-20)&
        # (odb_key_val_dic['col_2']>-75)&
        # (odb_key_val_dic['col_10']==1))
        
        
        ## Norht. Hem
        # zone='hn20'
        # ind_flag_tot = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']>20)&
        # (odb_key_val_dic['col_2']<75))
        
        # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']>20)&
        # (odb_key_val_dic['col_2']<75)&
        # (odb_key_val_dic['col_10']==1))
        
        # ##Glo
        zone='Glo'
        ind_flag_tot = np.where((odb_key_val_dic['col_1']<180)&
        (odb_key_val_dic['col_1']>-180))
        
        ind_flag = np.where((odb_key_val_dic['col_10']==1))
        

        
        #An depar
        data2plota =    -odb_key_val_dic['col_9'][ind_flag] ## col_7:obs SST, col_8 fg_depar-obs, col_9 an_depar-obs
        data2plota_tot = -odb_key_val_dic['col_9'][ind_flag_tot]      
        Nb=len(data2plota)
        Nb_tot=len(data2plota_tot)
        
        lon_bouees =    odb_key_val_dic['col_1'][ind_flag]
        lat_bouees =    odb_key_val_dic['col_2'][ind_flag]
        
        #Fg depar
        data2plotb =    -odb_key_val_dic['col_8'][ind_flag]
        data2plotb_tot = -odb_key_val_dic['col_8'][ind_flag_tot]  
        
        Nb_fg=len(data2plotb)
        Nb_tot_fg=len(data2plotb_tot)
        Nb_an=len(data2plota)

        #diurnal mean compute
        an_depar = data2plota.mean()
        fg_depar = data2plotb.mean()

        date_3 = date_dt.strftime("%d/%m/%y %H:%M")
        new_row = {'date':date_3, 'expe':expe[e],'andepar': an_depar,'fgdepar':fg_depar, 'Nbbouees':Nb, 'Nbbouees_tot':Nb_tot, 'Nbbouees_fg':Nb_fg, 'Nbbouees_an':Nb_an}
            
        # Use the loc method to add the new row to the DataFrame
        df_diurne_bouees.loc[len(df_diurne_bouees)] = new_row
  
    # drop first line
df_diurne_bouees.drop(index=df_diurne_bouees.index[0], axis=0, inplace=True) ## DELETE FIRST LINE
    
    ## Enregistre df separé
    # df_diurne_bouees.to_csv('DF_diurne'+'/'+'diurne_bouees_septembre_octobre_'+'.csv')
df_GL55_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[0]))]
df_GKOR_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[1]))]
df_GKPH_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[2]))]
df_GN84_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[3]))]
df_GO4A_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[4]))]
df_GOJQ_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[5]))]
df_GKCJ_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[6]))]


#%%
# FIGURE Eval SST, avec Mercator et bouees
# %matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np

fig1, ax = plt.subplots(1, 2, figsize=(28, 5.5),dpi=180)
fig1.tight_layout(pad=100)  
fig1.subplots_adjust(wspace=0.3)
L=len(df_GN84_bouees.date)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84_bouees.date[0:len(df_GN84_bouees.date):N]

# ax[0].set_ylabel('Analyse - Buoys (°C)',fontsize=23)
# ax[0].yaxis.set_tick_params(labelsize=23)
# ax[0].tick_params(labelbottom=False)
ax[0].grid(True, alpha=0.5)

## ARRANGE FIGURE
## BIAIS
ax[0].grid(True, alpha=0.5)
ax[0].set_xticks([])
ax[0].set_ylabel('Analyse - buoys (°C)',fontsize=24)
ax[0].yaxis.set_tick_params(labelsize=23)
xlabels = df_GN84_bouees.date[0:len(df_GN84_bouees.date):N]
ax[0].set_xticks(xind,labels=xlabels,rotation=45,fontsize=18)
ax[0].yaxis.set_tick_params(labelsize=24)

# plt.suptitle('-andepar,'+str(zone),fontsize=22, y=0.88)
ax[0].plot(df_GL55_bouees.date,df_GL55_bouees['andepar'],color=color_GL55,label='L.26 ',linewidth=2,linestyle='-')
ax[0].plot(df_GKOR_bouees.date,df_GKOR_bouees['andepar'],color=color_GKOR,label='L.12',linewidth=2,linestyle='-')
ax[0].plot(df_GKPH_bouees.date,df_GKPH_bouees['andepar'],color=color_GKPH,label='R.50 ',linewidth=2,linestyle='-')
ax[0].plot(df_GKCJ_bouees.date,df_GKCJ_bouees['andepar'],color=color_GKCJ,label='GKCJ',linewidth=2,linestyle='-')
ax[0].plot(df_GN84_bouees.date,df_GN84_bouees['andepar'],color=color_GN84,label='Xp.OML ',linewidth=2,linestyle='-')
# ax[0].plot(df_GO4A_bouees.date,df_GO4A_bouees['andepar'],color=color_GO4A,label='GO4A - buoys',linewidth=2,linestyle='-')
# ax[0].plot(df_GOJQ_bouees.date,df_GOJQ_bouees['andepar'],color=color_GOJQ,label='GPJQ - buoys',linewidth=2,linestyle='-')
ax[0].set_ylim(-0.2,0.35)
ax2 = ax[0].twinx()
ax2.plot(df_GL55_bouees.date,df_GL55_bouees['Nbbouees'],color=color_GL55,label='L.26 ',linewidth=1,linestyle=':')
ax2.plot(df_GKOR_bouees.date,df_GKOR_bouees['Nbbouees'],color=color_GKOR,label='L.12 ',linewidth=1,linestyle=':')
ax2.plot(df_GKPH_bouees.date,df_GKPH_bouees['Nbbouees'],color=color_GKPH,label='R.50',linewidth=1,linestyle=':')
ax2.plot(df_GKCJ_bouees.date, df_GKCJ_bouees['Nbbouees'], color=color_GKCJ, label='GKCJ ',linestyle=':', linewidth=1)
ax2.plot(df_GN84_bouees.date, df_GN84_bouees['Nbbouees'], color=color_GN84, label='Buoys numbers',linestyle=':', linewidth=1)
# ax2.plot(df_GO4A_bouees.date, df_GO4A_bouees['Nbbouees'], color=color_GO4A, label='Buoys numbers',linestyle=':', linewidth=1)
# ax2.plot(df_GOJQ_bouees.date, df_GOJQ_bouees['Nbbouees'], color=color_GOJQ, label='Buoys numbers',linestyle=':', linewidth=1)

ax2.set_ylabel('Buoys numbers',fontsize=23,color='steelblue')
ax2.yaxis.set_tick_params(labelsize=23,labelcolor='steelblue')
ax2.yaxis.set_label_coords(1.10, .5)
ax[0].yaxis.set_label_coords(-0.08, .5)

## ARRANGE FIGURE
## BIAIS
ax[1].grid(True, alpha=0.5)
ax[1].set_xticks([])
ax[1].set_ylabel('Forecast - buoys (°C)',fontsize=24)
ax[1].yaxis.set_tick_params(labelsize=23)
xlabels = df_GN84_bouees.date[0:len(df_GN84_bouees.date):N]
ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=18)
ax[1].yaxis.set_tick_params(labelsize=24)

# plt.suptitle('-fgdepar, '+str(zone),fontsize=22, y=0.88)
ax[1].plot(df_GL55_bouees.date,df_GL55_bouees['fgdepar'],color=color_GL55,label='L.26 ',linewidth=2,linestyle='-')
ax[1].plot(df_GKOR_bouees.date,df_GKOR_bouees['fgdepar'],color=color_GKOR,label='L.12 ',linewidth=2,linestyle='-')
ax[1].plot(df_GKPH_bouees.date,df_GKPH_bouees['fgdepar'],color=color_GKPH,label='R.50',linewidth=2,linestyle='-')
ax[1].plot(df_GKCJ_bouees.date,df_GKCJ_bouees['fgdepar'],color=color_GKCJ,label='Ref.noOML ',linewidth=2,linestyle='-')
ax[1].plot(df_GN84_bouees.date,df_GN84_bouees['fgdepar'],color=color_GN84,label='Xp.OML ',linewidth=2,linestyle='-')
# ax[1].plot(df_GO4A_bouees.date,df_GO4A_bouees['fgdepar'],color=color_GO4A,label='GO4A - buoys',linewidth=2,linestyle='-')
# ax[1].plot(df_GOJQ_bouees.date,df_GOJQ_bouees['fgdepar'],color=color_GOJQ,label='GOJQ - buoys',linewidth=2,linestyle='-')


ax2 = ax[1].twinx()
ax2.set_ylabel('Buoys numbers',fontsize=24,color='steelblue')
ax2.yaxis.set_tick_params(labelsize=24,labelcolor='steelblue')
ax2.axes.get_yaxis().set_visible(True)
ax2.plot(df_GL55_bouees.date, df_GL55_bouees['Nbbouees_fg'], color=color_GL55, label='L.26 ',linestyle=':', linewidth=1)
ax2.plot(df_GKOR_bouees.date, df_GKOR_bouees['Nbbouees_fg'], color=color_GKOR, label='L.12',linestyle=':', linewidth=1)
ax2.plot(df_GKPH_bouees.date, df_GKPH_bouees['Nbbouees_fg'], color=color_GKPH, label= 'Ref.noOML ',linestyle=':', linewidth=1)
ax2.plot(df_GN84_bouees.date, df_GN84_bouees['Nbbouees_fg'], color=color_GN84, label='Xp.OML',linestyle=':', linewidth=1)
# ax2.plot(df_GO4A_bouees.date, df_GO4A_bouees['Nbbouees_fg'], color=color_GO4A, label='GO4A - buoys',linestyle=':', linewidth=1)
# ax2.plot(df_GOJQ_bouees.date, df_GOJQ_bouees['Nbbouees_fg'], color=color_GOJQ, label='GOJQ - buoys',linestyle=':', linewidth=1)
ax2.yaxis.set_label_coords(1.11, .5) # place of labels
ax[1].yaxis.set_label_coords(-0.08, .5)
ax[1].set_ylim(-0.2,0.35)

# ax[0].axhline(y=0, color='k', linestyle='--', linewidth=1)
# ax[1].axhline(y=0, color='k', linestyle='--', linewidth=1)

ax[0].axhline(y=0, color='k', linestyle='-', linewidth=1, alpha=0.7)
ax[1].axhline(y=0, color='k', linestyle='-', linewidth=1, alpha=0.7)
ax[1].legend( loc='upper left',prop={'size': 21})

plt.show()
figname1 = dir_im+'an+fgdepar_firstexpe'+'_'+str(zone)+'_'+str(period)+'.png'
fig1.savefig(figname1,dpi=300, format='png',bbox_inches='tight')
# figname = dir_im+'fgdepar'+'_'+str(zone)+'.pdf'
# fig.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')
#%%
# # FIGURE Eval SST, avec Mercator et bouees
# %matplotlib inline
# import matplotlib.pyplot as plt
# plt.style.use('seaborn-white')
# import numpy as np


# fig2, ax = plt.subplots(1, 1, figsize=(16, 6),dpi=300)
# fig2.tight_layout(pad=8)  


# L=len(df_GN84_bouees.date)
# print(L)
# lst = list(np.arange(1,L+1))
# print(lst)
# xx=lst
# xind = xx[0:len(xx):N] # Pas des dates plot
# xlabels = df_GN84_bouees.date[0:len(df_GN84_bouees.date):N]

# ## ARRANGE FIGURE
# ## BIAIS
# ax.set_xticks([])
# ax.set_ylabel('Bias (°C)',fontsize=18)
# ax.yaxis.set_tick_params(labelsize=16)
# xlabels = df_GN84_bouees.date[0:len(df_GN84_bouees.date):N]
# ax.set_xticks(xind,labels=xlabels,rotation=45,fontsize=16)
# ax.set_xlabel('Date',fontsize=18)
# ax.yaxis.set_tick_params(labelsize=18)


# plt.suptitle('-fgdepar, '+str(zone),fontsize=22, y=0.88)

# ax.plot(df_GL55_bouees.date,df_GL55_bouees['fgdepar'],color=color_GL55,label='L.26 - buoys',linewidth=3,linestyle='--')
# ax.plot(df_GKOR_bouees.date,df_GKOR_bouees['fgdepar'],color=color_GKOR,label='L.12 - buoys',linewidth=3,linestyle='--')
# ax.plot(df_GKPH_bouees.date,df_GKPH_bouees['fgdepar'],color=color_GKPH,label='R.50 - buoys',linewidth=3,linestyle='--')
# ax.plot(df_GN84_bouees.date,df_GN84_bouees['fgdepar'],color=color_GN84,label='Xp.OML - buoys',linewidth=3,linestyle='--')
# ax.plot(df_GO4A_bouees.date,df_GO4A_bouees['fgdepar'],color=color_GO4A,label='GO4A - buoys',linewidth=3,linestyle='--')

# ax.plot(df_GKCJ_bouees.date,df_GKCJ_bouees['fgdepar'],color=color_GKCJ,label='GKCJ - buoys',linewidth=3,linestyle='--')
# ax.plot(df_GOJQ_bouees.date,df_GOJQ_bouees['fgdepar'],color=color_GOJQ,label='GOJQ - buoys',linewidth=3,linestyle='--')


# ax2 = ax.twinx()
# ax2.plot(df_GL55_bouees.date, df_GL55_bouees['Nbbouees_fg'], color=color_GL55, label='L.26 - buoys',linestyle=':', linewidth=2)
# ax2.plot(df_GKOR_bouees.date, df_GKOR_bouees['Nbbouees_fg'], color=color_GKOR, label='L.12 - buoys',linestyle=':', linewidth=2)
# ax2.plot(df_GKPH_bouees.date, df_GKPH_bouees['Nbbouees_fg'], color=color_GKPH, label='L.12 - buoys',linestyle=':', linewidth=2)
# ax2.plot(df_GN84_bouees.date, df_GN84_bouees['Nbbouees_fg'], color=color_GN84, label='Xp.OML - buoys',linestyle=':', linewidth=2)
# ax2.plot(df_GO4A_bouees.date, df_GO4A_bouees['Nbbouees_fg'], color=color_GO4A, label='GO4A - buoys',linestyle=':', linewidth=2)

# ax2.plot(df_GKCJ_bouees.date, df_GKCJ_bouees['Nbbouees_fg'], color=color_GKCJ, label='GKCJ - buoys',linestyle=':', linewidth=2)
# ax2.plot(df_GOJQ_bouees.date, df_GOJQ_bouees['Nbbouees_fg'], color=color_GOJQ, label='GOJQ - buoys',linestyle=':', linewidth=2)

# ax2.set_ylabel('Buoys numbers',fontsize=18,color='steelblue')
# ax2.yaxis.set_tick_params(labelsize=16,labelcolor='steelblue')

# ax.axhline(y=0, color='k', linestyle='-', linewidth=1)
# ax.grid(alpha=0.5)
# # ax.legend([ 'Xp.OMLe','SST.cst.mer','Xp.OMLe_bouees','SST.cst.bouées'],fontsize=18, loc='best')
# ax.legend(loc='best',prop={'size': 16})
# ax.legend(bbox_to_anchor=(1.00, 0), loc='upper left',prop={'size': 16})

# ax2.axes.get_yaxis().set_visible(True)
# plt.show()

# figname2 = dir_im+'fgdepar'+'_'+str(zone)+'_'+str(period)+'.png'
# fig2.savefig(figname2,dpi=300, format='png',bbox_inches='tight')
# # figname = dir_im+'fgdepar'+'_'+str(zone)+'.pdf'
# # fig.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')
