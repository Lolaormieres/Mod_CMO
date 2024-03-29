#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 17:59:55 2023

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
pathres='/home/ormieresl/Routines/DF_bouees_vf/'
dir_fig='/cnrm/recyf/Data/users/ormieresl/plot_evol_sst'

# experience to evaluate
# expe=['GN84','GO4A']
expe=['GKCJ']
version='double'

# reference
expe_mer=['GO48']
compar='prev-mer'+expe_mer[0]

# Levels numbers
nb_level=12

# Start date
date_init = datetime(2022, 9,2, 0)
deltat = timedelta(days=1)


# Days numbers
nb=1
nbDate=nb # buyos dates
period = timedelta(days=nb)

# échéances
ech_fin=102
ech_ini=0
ech_step=6

# cterm_for='96'
zone =  'nordat'

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

## Definition color
color_GL55='#648FFF'  ## GL55
color_GKOR='#785EF0'  ## GKOR
color_GKPH='red'  ## GKPH '#DC267F'
color_GM6E='#FE6100'  ## GM6E
color_GMOT='#009E73'  ## GMOT
color_GK7C='black' 
color_GKCJ='dimgrey' 
color_GN3C='gold'
color_GN84='red'
color_GNBA='lightseagreen'
color_GO4A='orange'
color_GNSR='lightseagreen'


########################################################
#definition des champs a extraire
########################################################
## ARP - assimilation
cfields = ''
cblock = 'surfan'
ckind = 'analysis'
cut = 'assim'
cmodel='surfex'
cterm='0'
formats='fa'
fill_cmo='surf'
# cblock='forecast'
# ckind='historic'

## ARP - production - prevision
cblock_for = 'forecast'
ckind_for = 'historic'
cmodel_for = 'surfex'
# cterm_for='96'
cterm_for='0'
cut_for = 'production'
formats_for='fa'
fill_for=''

## Mercator 
cblock_oc = 'c933'
ckind_oc = 'geofields'
cfields_oc = 'ocean'
#cfields_oc = 'sst'
cmodel_oc = 'arpege'
cut_oc = 'assim'

geometry = 'global1798'

# =============================================================================
# ## Lectur ficher fa
# =============================================================================
#LECTUR FILES
def get_file(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,
          'local':'tmpbouees_an_double.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpbouees_an_double.fic',openmode='r')
    os.remove('tmpbouees_an_double.fic')
    return r

#%%

ZONE='NORDAT'
Type='Analyse'

zone_GK7C = 'eurat'
zone='eurat'
zone_GN84 = 'eurat'
nom_zone = 'eurat'
# zone_GK7C = 'nord_atlantique'
# zone='nord-atlantique'
# zone_GN84 = 'nord-atlantique'
# nom_zone = 'Nord-Atlantique'
# zone_GK7C = 'tropiques'
# zone='tropiques'
# zone_GN84 = 'tropiques'
# nom_zone = 'tropiques'


#%%
## PLOT CYCLE DIURNE - BOUEES
expe=['GN84','GO4A']
XP='GN84'
SST_bouees=[]

### Ceate Dataframe
raw_data = {'date': [''],
            'lon': [''],
            'lat': [''],
        'expe': [''],
        'SST': [''],
        'Nbbouees':['']}
df_diurne_bouees= pd.DataFrame(raw_data, columns = ['date','lon','lat','expe','SST','Nbbouees'])

for e in range(len(expe)):
    ## BOUES BOUCLE TEST
    import datetime as dt
    year=2022
    month=9 # 4 pour Avril (et non 5 Mai)
    day=2
    hour=0
    Delta=6
    nbDate=106
    Date_debut = dt.datetime(year, month, day, hour)
    timedelta = dt.timedelta(0,0,0,0,0,Delta)
    toto=[Date_debut+i*timedelta for i in range(nbDate)]
    
    # #Répertoire de travail
    chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(expe[e])+'/'
    
    # BOUCLE lecture des fichiers
    # cterm_for= ['0','6','12','18','0','6','12','18','0','6','12','18','0','6','12','18','0']
    for date_dt in toto:
        print(date_dt)
        date = '%4.4i%2.2i%2.2iT%2.2i00A'%(date_dt.year,\
        date_dt.month,\
        date_dt.day,\
        date_dt.hour)
        print(date)
        # print(toto)
        date_2 = '%2.2i/%2.2i/%2.2i %2.2i'%(date_dt.year,\
        date_dt.month,\
        date_dt.day,\
        date_dt.hour)
        print(date_2)
        date_3 = '%2.2i%2.2i%2.2i%2.2i00'%(date_dt.year,\
        date_dt.month,\
        date_dt.day,\
        date_dt.hour)
        print(date_3)
        
        fic_nc = chemin_nc +'{}/surfan/odb-ecma.canari.surf.nc'.format(date)
        print(fic_nc)
        if not os.path.isfile(fic_nc):
            print('pas de fichier '+fic_nc)
        #print(len(fic_nc))   
    
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
        
        # EURAT
        ind_flag = np.where((odb_key_val_dic['col_1']>-35)&
        (odb_key_val_dic['col_1']<45)&
        (odb_key_val_dic['col_2']>20)&
        (odb_key_val_dic['col_2']<72))
        
        #N ATL
        # ind_flag = np.where((odb_key_val_dic['col_1']<0)&
        # (odb_key_val_dic['col_1']>-80)&
        # (odb_key_val_dic['col_2']>20)&
        # (odb_key_val_dic['col_2']<70))
        
        # TROPIQUES
        # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']>-20)&
        # (odb_key_val_dic['col_2']<20))
        
        ## Norht. Hem
        # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']>20)&
        # (odb_key_val_dic['col_2']<75))
        
              
        ## South. Hem
        # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
        # (odb_key_val_dic['col_1']>-180)&
        # (odb_key_val_dic['col_2']<-20)&
        # (odb_key_val_dic['col_2']>-75)&
        #  (odb_key_val_dic['col_2']==1))
        
        #data2plota = odb_key_val_dic['obsvalue'][ind_flag]
        data2plota =    odb_key_val_dic['col_9'][ind_flag] ## col_7:obs SST, col_8 fg_depar-obs, col_9 an_depar-obs
        lon_bouees =    odb_key_val_dic['col_1'][ind_flag]
        lat_bouees =    odb_key_val_dic['col_2'][ind_flag]
        datamin   = 275
        datamax   = 300
        
        x_bouees_pt = data2plota.mean()
        Nb=len(data2plota)
        print('data2plota',data2plota)
        print('len sst bouees', len(data2plota))
        print('len lon bouees', len(lon_bouees))
        print('len lat bouees', len(lat_bouees))
        
        print('sst bouees', x_bouees_pt)
        print('lon bouees', lon_bouees)
        print('lat bouees', lat_bouees)
        
    
        print('SST xpt',x_bouees_pt)
        SST_bouees=np.append(SST_bouees,x_bouees_pt)
        print('SST append bouees', SST_bouees)
        # print('DATE bouees',date)
    
        # x = datetime(2022, 9, 8, 0, 0)
        # print(x)
        # print(x.strftime("%d%m%y%H%M"))
        # x2=x.strftime("%d%m%y%H%M")
        
        # print(date_dt.strftime("%d%m%y%H%M"))
        # date_3 = date_dt.strftime("%d%m%y%H%M")
        
        date_3 = date_dt.strftime("%d/%m/%y %H:%M")
        # date_3 = date_dt.strftime("%d/%m/%y")
        # date_3=int(date_3)
        # print(type(date_3))
    
        new_row = {'date':date_3, 'lon':lon_bouees, 'lat':lat_bouees, 'expe':expe[e],'SST': x_bouees_pt, 'Nbbouees':Nb}
            
        # Use the loc method to add the new row to the DataFrame
        df_diurne_bouees.loc[len(df_diurne_bouees)] = new_row
    # print(df_diurne_bouees)
    # print(df_diurne_bouees.iloc[1])
    
    
    # drop first line
df_diurne_bouees.drop(index=df_diurne_bouees.index[0], axis=0, inplace=True) ## DELETE FIRST LINE
    
    ## Enregistre df separé
    # df_diurne_bouees.to_csv('DF_diurne'+'/'+'diurne_bouees_septembre_octobre_'+'.csv')
df_GN84_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[0]))]
df_GO4A_bouees=df_diurne_bouees.loc[(df_diurne_bouees['expe']==str(expe[1]))]

#%% OPER
## PLOT CYCLE DIURNE - BOUEES
expe=['GN84','GN3C']
XP='GKCJ'
SST_bouees=[]

### Ceate Dataframe
raw_data = {'date': [''],
            'lon': [''],
            'lat': [''],
        'expe': [''],
        'SST': [''],
        'Nbbouees':['']}
df_diurne_bouees_oper= pd.DataFrame(raw_data, columns = ['date','lon','lat','expe','SST','Nbbouees'])



## BOUES BOUCLE TEST
import datetime as dt
year=2022
month=9 # 4 pour Avril (et non 5 Mai)
day=2
hour=0
Delta=6
nbDate=106
Date_debut = dt.datetime(year, month, day, hour)

timedelta = dt.timedelta(0,0,0,0,0,Delta)

toto=[Date_debut+i*timedelta for i in range(nbDate)]

reseaux=""

# #Répertoire de travail
chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+str(XP)+'/'

# BOUCLE lecture des fichiers
# cterm_for= ['0','6','12','18','0','6','12','18','0','6','12','18','0','6','12','18','0']
for date_dt in toto:
    print(date_dt)
    date = '%4.4i%2.2i%2.2iT%2.2i00A'%(date_dt.year,\
    date_dt.month,\
    date_dt.day,\
    date_dt.hour)
    print(date)
    # print(toto)
    date_2 = '%2.2i/%2.2i/%2.2i %2.2i'%(date_dt.year,\
    date_dt.month,\
    date_dt.day,\
    date_dt.hour)
    print(date_2)
    date_3 = '%2.2i%2.2i%2.2i%2.2i00'%(date_dt.year,\
    date_dt.month,\
    date_dt.day,\
    date_dt.hour)
    print(date_3)
    

    fic_nc = chemin_nc +'{}/surfan/odb-ecma.canari.surf.nc'.format(date)
    print(fic_nc)
    if not os.path.isfile(fic_nc):
        print('pas de fichier '+fic_nc)
    #print(len(fic_nc))   

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
    
    # EURAT
    ind_flag = np.where((odb_key_val_dic['col_1']>-35)&
    (odb_key_val_dic['col_1']<45)&
    (odb_key_val_dic['col_2']>20)&
    (odb_key_val_dic['col_2']<72))
    
    ##N ATL
    # ind_flag = np.where((odb_key_val_dic['col_1']<0)&
    # (odb_key_val_dic['col_1']>-80)&
    # (odb_key_val_dic['col_2']>20)&
    # (odb_key_val_dic['col_2']<70))
    
    ## TROPIQUES
    # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
    # (odb_key_val_dic['col_1']>-180)&
    # (odb_key_val_dic['col_2']>-20)&
    # (odb_key_val_dic['col_2']<20))
    
    ## Hem. Nord
    # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
    # (odb_key_val_dic['col_1']>-180)&
    # (odb_key_val_dic['col_2']>20)&
    # (odb_key_val_dic['col_2']<75))
    
                  
    ## South. Hem
    # ind_flag = np.where((odb_key_val_dic['col_1']<180)&
    # (odb_key_val_dic['col_1']>-180)&
    # (odb_key_val_dic['col_2']<-20)&
    # (odb_key_val_dic['col_2']>-75)&
    # (odb_key_val_dic['col_10']==1))
    
    #data2plota = odb_key_val_dic['obsvalue'][ind_flag]
    data2plota =    odb_key_val_dic['col_9'][ind_flag] ## col_7:obs SST, col_8 fg_depar-obs, col_9 an_depar-obs
    lon_bouees =    odb_key_val_dic['col_1'][ind_flag]
    lat_bouees =    odb_key_val_dic['col_2'][ind_flag]
    
    datamin   = 275
    datamax   = 300
    
    x_bouees_pt = data2plota.mean()
    Nb=len(data2plota)
    print('data2plota',data2plota)
    print('len sst bouees', len(data2plota))
    print('len lon bouees', len(lon_bouees))
    print('len lat bouees', len(lat_bouees))
    
    print('sst bouees', x_bouees_pt)
    print('lon bouees', lon_bouees)
    print('lat bouees', lat_bouees)
    

    print('SST xpt',x_bouees_pt)
    SST_bouees=np.append(SST_bouees,x_bouees_pt)
    print('SST append bouees', SST_bouees)
    # print('DATE bouees',date)

    # x = datetime(2022, 9, 8, 0, 0)
    # print(x)
    # print(x.strftime("%d%m%y%H%M"))
    # x2=x.strftime("%d%m%y%H%M")
    
    # print(date_dt.strftime("%d%m%y%H%M"))
    # date_3 = date_dt.strftime("%d%m%y%H%M")
    
    date_3 = date_dt.strftime("%d/%m/%y %H:%M")
    # date_3 = date_dt.strftime("%d/%m/%y")
    # date_3=int(date_3)
    # print(type(date_3))

    new_row = {'date':date_3, 'lon':lon_bouees, 'lat':lat_bouees, 'expe':XP,'SST': x_bouees_pt, 'Nbbouees':Nb}
        
    # Use the loc method to add the new row to the DataFrame
    df_diurne_bouees_oper.loc[len(df_diurne_bouees_oper)] = new_row
# print(df_diurne_bouees)
# print(df_diurne_bouees.iloc[1])


# drop first line
df_diurne_bouees_oper.drop(index=df_diurne_bouees_oper.index[0], axis=0, inplace=True) ## DELETE FIRST LINE

# ## Enregistre df separé
# df_diurne_bouees.to_csv('DF_diurne_oper'+'/'+'diurne_bouees_septembre_octobre_'+'.csv')


#%%
# FIGURE Eval SST, avec Mercator et bouees
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
# N=2

fig2, ax = plt.subplots(1, 1, figsize=(16, 6),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig2.tight_layout(pad=8)  

# Plot dates axe x
N=5

L=len(df_GN84.date)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84.date[0:len(df_GN84.date):N]

## ARRANGE FIGURE
## BIAIS
ax.set_xticks([])
ax.set_ylabel('Bias (°C)',fontsize=18)
ax.yaxis.set_tick_params(labelsize=16)
xlabels = df_GN84.date[0:len(df_GN84.date):N]
ax.set_xticks(xind,labels=xlabels,rotation=45,fontsize=16)
ax.set_xlabel('Date',fontsize=18)
ax.yaxis.set_tick_params(labelsize=18)

# plt.suptitle('SST ARPEGE analysée - SST Mercator (REF) \n\n',fontsize=18)
# plt.suptitle('SST analyse modèle CMO 1D - SST Mercator (REF), \n '+str(nom_zone),fontsize=18, y=0.95)
# plt.suptitle('Evaluation of the analysed SST of the OML model (Xp.final), \n zone: '+str(nom_zone),fontsize=22, y=0.95)
plt.suptitle('Evaluation of the analysed SST of the CMO model (Xp.final)',fontsize=22, y=0.97)
# ax[0].plot(df_GL55.date,df_GL55.biais,label='L.26',color=color_GL55,linewidth=3)
# ax[0].plot(df_GKOR.date,df_GKOR.biais,label='L.12',color=color_GKOR,linewidth=3)
# ax[0].plot(df_GKPH.date,df_GKPH.biais,label='R.50',color=color_GKPH,linewidth=3)
# ax[0].plot(df_GM6E.date,df_GM6E.biais,label='no.curr',color=color_GM6E,linewidth=3)
# ax[0].plot(df_GMOT.date,df_GMOT.biais,label='no.curr.bathy',color=color_GMOT,linewidth=3)
ax.plot(df_GN84.date,df_GN84.biais,label='Xp.finale - Mercator',color=color_GN84,linewidth=3)
# ax[0].plot(df_OPER.date,df_OPER.biais,label='OPER.OSTIA',color=color_OPER)
ax.plot(df_GN84.date,df_GKCJ.biais,label='oper.2024 - Mercator',color=color_GKCJ,linewidth=3)
# ax.plot(df['date'],-df['SST'],color=color_GN84,label='GN84',linewidth=1,linestyle=':')
ax.plot(df_GN84.date,-df_GN84_bouees['SST'],color=color_GN84,label='Xp.final - buoys',linewidth=3,linestyle='--')
ax.plot(df_GO4A_bouees.date,-df_GO4A_bouees['SST'],color=color_GO4A,label='GO4A - buoys',linewidth=3,linestyle='--')
ax.plot(df_GN84.date,-df_diurne_bouees_oper['SST'],label='oper.2024 - buoys',color=color_GKCJ,linewidth=3,linestyle='--')



ax2 = ax.twinx()
ax2.plot(df_GN84_bouees.date, df_GN84_bouees['Nbbouees'], color='steelblue', label='Buoys numbers',linestyle=':', linewidth=2)
# ax2.plot(df_GN84.date, df_diurne_bouees_oper['Nbbouees'], color='aqua', label='Buoys numbers',linestyle=':', linewidth=2)
ax2.set_ylabel('Buoys numbers',fontsize=18,color='steelblue')
ax2.yaxis.set_tick_params(labelsize=16,labelcolor='steelblue')

ax.axhline(y=0, color='k', linestyle='-', linewidth=1)
ax.legend([ 'Xp.finale','SST.cst.mer','XP.finale_bouees','SST.cst.bouées'],fontsize=18, loc='best')
# ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)
ax.legend(loc='best',prop={'size': 16})
ax.legend(bbox_to_anchor=(1.00, 0), loc='upper left',prop={'size': 16})



#%%
## STD
# xlabels = df_GL55.date[0:len(df_GL55.date):N]
# ax[1].set_xticks(xind,labels=xlabels,rotation=45,fontsize=16)
# ax[1].set_xlabel('Date',fontsize=16)
# ax[1].set_ylabel('Std (°C)',fontsize=16)
# ax[1].yaxis.set_tick_params(labelsize=16)
# # ax[1].plot(df_GKOR.date,df_GKOR.eqm,label='GKOR',color=color_GKOR,linewidth=3)
# # ax[1].plot(df_GL55.date,df_GL55.eqm,label='GL55',color=color_GL55,linewidth=3)
# # ax[1].plot(df_GKPH.date,df_GKPH.eqm,label='GKPH',color=color_GKPH,linewidth=3)
# # ax[1].plot(df_GM6E.date,df_GM6E.eqm,label='GM6E',color=color_GM6E,linewidth=3)
# # ax[1].plot(df_GMOT.date,df_GMOT.eqm,label='GMOT',color=color_GMOT,linewidth=3)
# ax[1].plot(df_GN84.date,df_GN84.eqm,label='Xp.finale',color=color_GN84,linewidth=3)
# # ax[1].plot(df_OPER.date,df_OPER.eqm,label='OPER.OSTIA',color=color_OPER)
# ax[1].plot(df_GK7C.date,df_GK7C.eqm,label='MERCATOR moyen',color=color_GK7C,linewidth=3)
# ax[1].legend(['L.12','L.26','R.50','no.curr','no.curr.bathy', 'Xp.finale','SST.cst.mer'],fontsize=14, loc='upper right')
# ax[1].legend([ 'Xp.finale','SST.cst.mer','XP.finale_bouees'],fontsize=14, loc='upper right')
# plt.show()


# figname = '/cnrm/recyf/Data/users/ormieresl/SST_evol_temporelle_analyse-mercator_'+'_'+str(zone)+'.png'
# fig2.savefig(figname,dpi=300, format='png',bbox_inches='tight')
# figname = '/cnrm/recyf/Data/users/ormieresl/SST_evol_temporelle_analyse-mercator_'+'_'+str(zone)+'.pdf'
# fig2.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')


#%%
# FIGURE Eval SST, avec Mercator et bouees
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
# N=2

fig2, ax = plt.subplots(1, 1, figsize=(16, 6),dpi=300)
# fig, ax = plt.subplots(2, 1)
fig2.tight_layout(pad=8)  

# Plot dates axe x
N=5

L=len(df_GN84_bouees.date)
print(L)
lst = list(np.arange(1,L+1))
print(lst)
xx=lst
xind = xx[0:len(xx):N] # Pas des dates plot
xlabels = df_GN84_bouees.date[0:len(df_GN84.date):N]

## ARRANGE FIGURE
## BIAIS
ax.set_xticks([])
ax.set_ylabel('Bias (°C)',fontsize=18)
ax.yaxis.set_tick_params(labelsize=16)
xlabels = df_GN84_bouees.date[0:len(df_GN84_bouees.date):N]
ax.set_xticks(xind,labels=xlabels,rotation=45,fontsize=16)
ax.set_xlabel('Date',fontsize=18)
ax.yaxis.set_tick_params(labelsize=18)

# plt.suptitle('SST ARPEGE analysée - SST Mercator (REF) \n\n',fontsize=18)
# plt.suptitle('SST analyse modèle CMO 1D - SST Mercator (REF), \n '+str(nom_zone),fontsize=18, y=0.95)
# plt.suptitle('Evaluation of the analysed SST of the OML model (Xp.final), \n zone: '+str(nom_zone),fontsize=22, y=0.95)
plt.suptitle('Evaluation of the analysed SST of the CMO model (Xp.final), zone: '+str(ZONE)+' '+str(Type),fontsize=22, y=0.97)
# ax[0].plot(df_GL55.date,df_GL55.biais,label='L.26',color=color_GL55,linewidth=3)
# ax[0].plot(df_GKOR.date,df_GKOR.biais,label='L.12',color=color_GKOR,linewidth=3)
# ax[0].plot(df_GKPH.date,df_GKPH.biais,label='R.50',color=color_GKPH,linewidth=3)
# ax[0].plot(df_GM6E.date,df_GM6E.biais,label='no.curr',color=color_GM6E,linewidth=3)
# ax[0].plot(df_GMOT.date,df_GMOT.biais,label='no.curr.bathy',color=color_GMOT,linewidth=3)
# ax.plot(df_GN84.date,df_GN84.biais,label='Xp.finale - Mercator',color=color_GN84,linewidth=3)
# ax[0].plot(df_OPER.date,df_OPER.biais,label='OPER.OSTIA',color=color_OPER)
# ax.plot(df_GN84.date,df_GKCJ.biais,label='oper.2024 - Mercator',color=color_GKCJ,linewidth=3)
# ax.plot(df['date'],-df['SST'],color=color_GN84,label='GN84',linewidth=1,linestyle=':')
ax.plot(df_GN84_bouees.date,-df_GN84_bouees['SST'],color=color_GN84,label='Xp.final - buoys',linewidth=3,linestyle='--')
ax.plot(df_GO4A_bouees.date,-df_GO4A_bouees['SST'],color=color_GO4A,label='GO4A - buoys',linewidth=3,linestyle='--')
ax.plot(df_diurne_bouees_oper.date,-df_diurne_bouees_oper['SST'],label='oper.2024 - buoys',color=color_GKCJ,linewidth=3,linestyle='--')



ax2 = ax.twinx()
ax2.plot(df_GN84_bouees.date, df_GO4A_bouees['Nbbouees'], color='steelblue', label='Buoys numbers',linestyle=':', linewidth=2)
ax2.plot(df_GN84_bouees.date, df_diurne_bouees_oper['Nbbouees'], color='dimgrey', label='Buoys numbers',linestyle=':', linewidth=2)

# ax2.plot(df_GN84.date, df_diurne_bouees_oper['Nbbouees'], color='aqua', label='Buoys numbers',linestyle=':', linewidth=2)
ax2.set_ylabel('Buoys numbers',fontsize=18,color='steelblue')
ax2.yaxis.set_tick_params(labelsize=16,labelcolor='steelblue')

ax.axhline(y=0, color='k', linestyle='-', linewidth=1)
ax.legend([ 'Xp.finale','SST.cst.mer','XP.finale_bouees','SST.cst.bouées'],fontsize=18, loc='best')
# ax[0].legend(['L.26','L.12','R.50','no.curr','no.curr.bathy', 'SST.cst.mer'],fontsize=12)
ax.legend(loc='best',prop={'size': 16})
ax.legend(bbox_to_anchor=(1.00, 0), loc='upper left',prop={'size': 16})
