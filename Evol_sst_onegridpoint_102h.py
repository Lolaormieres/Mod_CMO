#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:31:25 2023
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


##CMO
name=['SFX.SST']
name_sal=['SFX.SALT_OC1']#,'SFX.SALT_OC2','SFX.SALT_OC3','SFX.SALT_OC4','SFX.SALT_OC5','SFX.SALT_OC6','SFX.SALT_OC7','SFX.SALT_OC8','SFX.SALT_OC9','SFX.SALT_OC10','SFX.SALT_OC11','SFX.SALT_OC12']#,'SFX.SALT_OC13','SFX.SALT_OC14','SFX.SALT_OC15','SFX.SALT_OC16','SFX.SALT_OC17','SFX.SALT_OC18','SFX.SALT_OC19','SFX.SALT_OC20''SFX.SALT_OC21','SFX.SALT_OC22','SFX.SALT_OC23','SFX.SALT_OC24','SFX.SALT_OC25','SFX.SALT_OC26']
name_u=['SFX.UCUR_OC1','SFX.UCUR_OC2','SFX.UCUR_OC3','SFX.UCUR_OC4','SFX.UCUR_OC5','SFX.UCUR_OC6','SFX.UCUR_OC7','SFX.UCUR_OC8','SFX.UCUR_OC9','SFX.UCUR_OC10','SFX.UCUR_OC11','SFX.UCUR_OC12','SFX.UCUR_OC13','SFX.UCUR_OC14','SFX.UCUR_OC15','SFX.UCUR_OC16','SFX.UCUR_OC17','SFX.UCUR_OC18','SFX.UCUR_OC19','SFX.UCUR_OC20','SFX.UCUR_OC21','SFX.UCUR_OC22','SFX.UCUR_OC23','SFX.UCUR_OC24','SFX.UCUR_OC25','SFX.UCUR_OC26']
name_v=['SFX.VCUR_OC1','SFX.VCUR_OC2','SFX.VCUR_OC3','SFX.VCUR_OC4','SFX.VCUR_OC5','SFX.VCUR_OC6','SFX.VCUR_OC7','SFX.VCUR_OC8','SFX.VCUR_OC9','SFX.VCUR_OC10','SFX.VCUR_OC11','SFX.VCUR_OC12','SFX.VCUR_OC13','SFX.VCUR_OC14','SFX.VCUR_OC15','SFX.VCUR_OC16','SFX.VCUR_OC17','SFX.VCUR_OC18','SFX.VCUR_OC19','SFX.VCUR_OC20','SFX.VCUR_OC21','SFX.VCUR_OC22','SFX.VCUR_OC23','SFX.VCUR_OC24','SFX.VCUR_OC25','SFX.VCUR_OC26']
## MERCATOR
name_mer=['SURF.THETAO1']#],'SURF.THETAO12']#,'SURF.THETAO2','SURF.THETAO3','SURF.THETAO4','SURF.THETAO5','SURF.THETAO6','SURF.THETAO7','SURF.THETAO8','SURF.THETAO9','SURF.THETAO10','SURF.THETAO11','SURF.THETAO12','SURF.THETAO13','SURF.THETAO14','SURF.THETAO15','SURF.THETAO16','SURF.THETAO17','SURF.THETAO18','SURF.THETAO19','SURF.THETAO20','SURF.THETAO21','SURF.THETAO22','SURF.THETAO23','SURF.THETAO24','SURF.THETAO25','SURF.THETAO26']
name_mer_sal=['SURF.SALINO1','SURF.SALINO2','SURF.SALINO3','SURF.SALINO4','SURF.SALINO5','SURF.SALINO6','SURF.SALINO7','SURF.SALINO8','SURF.SALINO9','SURF.SALINO10','SURF.SALINO11','SURF.SALINO12']#,'SURF.SALINO13','SURF.SALINO14','SURF.SALINO15','SURF.SALINO16','SURF.SALINO17','SURF.SALINO18','SURF.SALINO19','SURF.SALINO20','SURF.SALINO21','SURF.SALINO22','SURF.SALINO23','SURF.SALINO24','SURF.SALINO25','SURF.SALINO26']
## ATMO
name_atm=[{'shortName':'ssrd'},{'shortName':'10u'},{'shortName':'10v'}] #,{'shortNameECMF':'100v'}] 
# name_atm=[{'shortName':'10u'}] 
# name_atm=[{'shortName':'10v'}] 

# Lecture des niveaux de la CMO
nb_level=12

## Lecture des niveaux de la CMO
filename = 'level2.txt'
levelfic = np.loadtxt(filename, delimiter=',', dtype=float)
print('level', levelfic)
print(type(levelfic))
print(np.shape(levelfic))
level = levelfic[0:nb_level]

expe='GL55'
expe2='GL55'
expe3='GKPH'
xpref='GKOR'
#geometry='global798c22'
geometry='global1798'
#chemin='/d0/images/ormieresl/'
chemin='/home/ormieresl/Scripts/'


Dir_im='/cnrm/recyf/Data/users/ormieresl/plot_sst_forecast_102h_testcoeff/'

## ARP
cfields=''
cblock='surfan'
ckind='analysis'
cut='assim'
formats='fa'
fill_cmo='surf'
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

cterm_6='6'
cterm_12='12'
cterm_18='18'
cterm_24='24'
cterm_48='48'
cterm_54='54'
cterm_60='60'
cterm_p6='6'

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

## ARP Atmo
fill_atm=''
grille="glob01"
geometry_atm=grille
cblock_atm='forecast'
ckind_atm='historic'
nbmb=1
formats_atm='grib'
ckind_atm='gridpoint'
cblock_atm='forecast'
cfields_atm='grid'
cmodel_atm='arpege'
cut_atm='production'

## Definition color Pale
color_GL55='steelblue' ## GL55 '#648FFF'
color_GKOR='hotpink'  ## GKOR
color_GKPH= 'mediumpurple' ## GKPH '#DC267F'
color_GM6E='paleturquoise'
color_GMOT='peachpuff'
color_GN84='palegreen'#'r'
color_GO4A='rosybrown'
color_GNSR='palegreen'
color_GNYG='yellow'
color_GOJQ='orange'
color_GOJQ='peachpuff'
color_GOJQ='midnightblue'
color_GOJQ='rosybrown'
color_GO4A='steelblue'

color_expe=color_GMOT
color_mer='gold'

##Lat/lon du pt de grille a etudier
# long=5
# lati=41

# long=7
# lati=-31

# long=4
# lati=38
long=3
lati=37

# long=-9
# lati=18

# long=5
# lati=38

# long=7
# lati=42

# long=-22
# lati=10

# lonmin=0
# lonmax=21

# latmin=30
# latmax=50


#%%read fa 
def get_file(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,
          'local':'tmp102.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmp102.fic',openmode='r')
    os.remove('tmp102.fic')
    return r

def get_file_atm(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,
          'local':'tmpatm.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr', 'origin':'historic',
          'nativefmt':'grib','now':True}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpatm.fic',openmode='r')
    return r
   
#%%
##sst Mercator
expe=['GO48']
nb=17
date_init = datetime(2022, 8,2 , 18)
date = date_init
dt = timedelta(hours=6)
date_valid = date+dt

##create dataframe
raw_data = {'date': [''],
            'lon': [''],
            'lat': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'ech': [''],
                'SST': ['']}
df_diurne_oc = pd.DataFrame(raw_data, columns = ['date','lon','lat','expe','var','level','ech','SST'])

for t in range(nb):
    date = date + dt
    print(date)
    for e in range(len(expe)):
        print(expe[e])        
        Pas='6h'

        cterm_oc= ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']#,  '6', '12', '18']# '24', '30', '36', '42', '48', '54', '60', '66', '72', '78', '84', '90', '96']
        cterm_oc= ['0']
        l=np.arange(0,102,6)
        c_trm_for=str(l)
        cterm_for=['0',  '6', '12', '18','24', '30', '36', '42', '48', '54', '60', '66', '72', '78', '84', '90', '96']
        for h in range(len(cterm_oc)):
            print(cterm_oc[h])
            try:
                #r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm)
                #r_ref=get_file(xpref,date,cblock,ckind,cfields,cmodel,cut,cterm)
                #r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for[h])
                #r_drown=get_file(expe_drown,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm)
                r_mer=get_file(expe,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc[h],geometry,formats,fill_oc)
                for i in range(len(name_mer)):
                    print(i,name_mer[i])
                    Y=r_mer.geometry.get_lonlat_grid()  #lon et lat
                    lon_mer=Y[0]
                    lat_mer=Y[1]

                    mask=(lon_mer<long)&(lon_mer>long-1)&(lat_mer>lati)&(lat_mer<lati+1) #masques sur lon et lat
                    
                    x=r_mer.readfield(name_mer[i])
                    crs=ccrs.PlateCarree()
                    
                    #1 zone
                    x_total=x.data[mask]
                    
                    #1 pt
                    x_pt=x_total[0]-273.15
                    print('SST en un point',x_pt)
                    print('lon MER',np.round(lon_mer[mask][0],2))
                    print('lat MER',np.round(lat_mer[mask][0],2))
                    
                    x_pt=np.array(x_pt)
      
                    date_2 = date.strftime("%d/%m/%y %H:%M")
    
                    new_row = {'date':date_2,'lon':lon_mer[0], 'lat':lat_mer[0], 'expe':expe[e], 'var':name_mer[i],'level':level[i], 'ech':cterm_for[h],'SST': x_pt}    
                      
                    # Use the loc method to add the new row to the DataFrame
                    df_diurne_oc.loc[len(df_diurne_oc)] = new_row
            except:
                    print('Pas de fichier')

# SUPRIMER 1er COL DF
df_diurne_oc.drop(index=df_diurne_oc.index[0], axis=0, inplace=True) ## DELETE FIRST LINE                         

#%%
## SST CMO analyse
expe=['GO4A']
nb=17
date_init = datetime(2022, 8,2, 18)
date = date_init
dt = timedelta(hours=6)
date_valid = date+dt

##create dataframe
raw_data = {'date': [''],
            'lon': [''],
            'lat': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'ech': [''],
                'SST': ['']}
df_diurne_an = pd.DataFrame(raw_data, columns = ['date','lon','lat','expe','var','level','ech','SST'])

for t in range(nb):
    date = date + dt
    print(date)
    for e in range(len(expe)):
        print(expe[e])
        Pas='6h'
        cterm_oc= ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']#,  '6', '12', '18']# '24', '30', '36', '42', '48', '54', '60', '66', '72', '78', '84', '90', '96']
        cterm_oc= ['0']

        cterm_for=['0',  '6', '12', '18','24', '30', '36', '42', '48', '54', '60', '66', '72', '78', '84', '90', '96']
        for h in range(len(cterm_oc)):
            print(cterm_oc[h])
            # try:
            #r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm)
            #r_ref=get_file(xpref,date,cblock,ckind,cfields,cmodel,cut,cterm)
            #r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for[h])

            # r_mer=get_file(expe,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc[h],geometry,formats,fill_oc)
            r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm_oc[h],geometry,formats,fill_cmo)
            # r_mer=get_file(expe,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc[h],geometry,formats)
            #SST=np.array(SST)
             
            for i in range(len(name)):
                print(i,name[i])
                
                Y=r.geometry.get_lonlat_grid()  #lon et lat
                lon=Y[0]
                lat=Y[1]

                mask=(lon<long)&(lon>long-1)&(lat>lati)&(lat<lati+1) #masques sur lon et lat
                x=r.readfield(name[i])
                crs=ccrs.PlateCarree()
                
                #1 zine
                x_total=x.data[mask]
                
                #1 pt
                x_pt=x_total[0]-273.15
                print('SST en un point',x_pt)
                print('lon expe',np.round(lon[mask][0],2))
                print('lat exp',np.round(lat[mask][0],2))
                
                x_pt=np.array(x_pt)
  
                date_2 = date.strftime("%d/%m/%y %H:%M")
                new_row = {'date':date_2,'lon':lon[0], 'lat':lat[0], 'expe':expe[e], 'var':name[i],'level':level[i], 'ech':cterm_for[h],'SST': x_pt}    
                # Use the loc method to add the new row to the DataFrame
                df_diurne_an.loc[len(df_diurne_an)] = new_row
            # except:
            #         print('Pas de fichier')

# SUPRIMER 1er COL DF
df_diurne_an.drop(index=df_diurne_an.index[0], axis=0, inplace=True) ## DELETE FIR

#%%
##sst CMO previsions
nb = 1
date_init = datetime(2022, 8, 2 , 0)
date = date_init
dt = timedelta(days=1)
date_valid = date+dt
period = timedelta(days=nb)
date_end = date+period    

##create dataframe
raw_data = {'date': [''],
            'lon': [''],
            'lat': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'ech': [''],
                'SST': ['']}
df_diurne = pd.DataFrame(raw_data, columns = ['date','lon','lat','expe','var','level','ech','SST'])

# expe=['GK7C','GL55','GKOR','GKPH','GM6E','GMOT', 'GN84','GNBA']
# expe=['GKCJ','GN84','GO4A','GOJQ','GGPN']


expe=['GKCJ','GN84','GPPS','GPPU','GPPN','GNBA']
for e in range(len(expe)):
    print(expe[e]) 
    Pas='6h'
    cterm_for= ['0',  '6', '12', '18', '24', '30', '36', '42', '48', '54', '60', '66', '72', '78', '84', '90', '96']
    for h in range(len(cterm_for)):
        print(cterm_for[h])
        try:
            #r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm)
            #r_ref=get_file(xpref,date,cblock,ckind,cfields,cmodel,cut,cterm)
            r_for=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for[h],geometry,formats,fill_cmo)
            #r_drown=get_file(expe_drown,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm)
            #r_mer=get_file(expe,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm)
            for i in range(len(name)):
                print(i,name[i])
                X=r_for.geometry.get_lonlat_grid()  #lon et lat
                lon=X[0]
                lat=X[1]
                        
                mask=(lon<long)&(lon>long-1)&(lat>lati)&(lat<lati+1) #masques sur lon et lat
                x=r_for.readfield(name[i])
                crs=ccrs.PlateCarree()
                #1 zone
                x_total=x.data[mask]
                # 1pt
                x_pt=x_total[0]-273.15
                lon_expe=np.round(lon[mask][0],2)
                lat_expe=np.round(lat[mask][0],2)
                print('lon expe',lon_expe)
                print('lat expe',lat_expe)
                print('SST en un point',x_pt)
               
                x_pt=np.array(x_pt)
                date_2 = date.strftime("%d/%m/%y")
                new_row = {'date':date_2,'lon':lon[0], 'lat':lat[0], 'expe':expe[e], 'var':name[i],'level':level[i], 'ech':cterm_for[h],'SST': x_pt}
                    
                # Use the loc method to add the new row to the DataFrame
                df_diurne.loc[len(df_diurne)] = new_row
        except:
                print('Pas de fichier')
 # SUPRIMER 1er COL DF
df_diurne.drop(index=df_diurne.index[0], axis=0, inplace=True) ## DELETE FIRST LINE               

#%%
##sst observee
expe=['GN84','GO48']
##create dataframe
raw_data = {'date': [''],
            'lon': [''],
            'lat': [''],
        'expe': [''],
        'SST': ['']}
df_diurne_bouees= pd.DataFrame(raw_data, columns = ['date','lon','lat','expe','SST'])

#boucle sur les dates pour obtenir le nom de repertoir dans lesquels le fichier Netcdf est
import datetime as dt
year=2022
month=9 # 4 pour Avril (et non 5 Mai)
day=2
hour=0
Delta=6
nbDate=17
Date_debut = dt.datetime(year, month, day, hour)
timedelta = dt.timedelta(0,0,0,0,0,Delta)
toto=[Date_debut+i*timedelta for i in range(nbDate)]
reseaux=""

#répertoire
chemin_nc = '/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/GN84/'
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
    
    #mask lon/lat + anflag (col 10)
    ind_flag = np.where((odb_key_val_dic['col_1']>long-2)&
    (odb_key_val_dic['col_1']<long+2)&
    (odb_key_val_dic['col_2']>lati-2)&
    (odb_key_val_dic['col_2']<lati+2)&
    (odb_key_val_dic['col_10']==1))
        
    #data2plota = odb_key_val_dic['obsvalue'][ind_flag]
    data2plota =    odb_key_val_dic['col_7'][ind_flag]
    lon_bouees =    odb_key_val_dic['col_1'][ind_flag]
    lat_bouees =    odb_key_val_dic['col_2'][ind_flag]
    anflag_val =    odb_key_val_dic['col_10']
    
    datamin   = 275
    datamax   = 300
    
    x_bouees_pt = data2plota[0] - 273.15
    print('len sst bouees', len(data2plota))

    print('sst bouees', x_bouees_pt)
    lon_bouees_pt=np.round(lon_bouees[0],2)
    lat_bouees_pt=np.round(lat_bouees[0],2)
    print('lon bouees', lon_bouees_pt)
    print('lat bouees', lat_bouees_pt)

    date_3 = date_dt.strftime("%d%m%y%H%M")
    new_row = {'date':date_3, 'lon':lon_bouees_pt, 'lat':lat_bouees_pt, 'expe':expe[0],'SST': x_bouees_pt}
    # Use the loc method to add the new row to the DataFrame
    df_diurne_bouees.loc[len(df_diurne_bouees)] = new_row

# SUPRIMER 1er COL DF
df_diurne_bouees.drop(index=df_diurne_bouees.index[0], axis=0, inplace=True) ## DELETE FIRST LINE

# #%%
# SST_GKCJ = df_diurne.loc[((df_diurne['expe']=='GKCJ') )][['date','expe','var','level','ech', 'SST']]
# GKCJ=SST_GKCJ.SST.iloc[1]
# GKCJ=np.linspace(GKCJ,GKCJ,17)

# SST_GN84 = df_diurne.loc[((df_diurne['expe']=='GN84'))][['date','expe','var','level','ech', 'SST']] 
# #SST_GL55 = df_diurne.loc[((df_diurne['expe']=='GL55'))][['date','expe','var','level','ech', 'SST']] 
# SST_GO4A = df_diurne.loc[((df_diurne['expe']=='GO4A'))][['date','expe','var','level','ech', 'SST']] 
# SST_GOJQ = df_diurne.loc[((df_diurne['expe']=='GOJQ'))][['date','expe','var','level','ech', 'SST']] 

# SST_GO4A_an = df_diurne_an.loc[((df_diurne_an['expe']=='GO4A'))][['date','expe','var','level','ech', 'SST']]

# fig, ax = plt.subplots(nrows = 1,ncols = 1, sharex = False, figsize=(8, 5),dpi=150) 

# ## BIAIS
# # ax.set_xticks([])
# ax.set_ylabel('SST (°C)',fontsize=19)
# # xlabels = df_diurne_atm.ech[0:len(df_diurne_atm.ech):N]
# # ax.set_xticks(xind,labels=xlabels,fontsize=15)

# ax.set_xlabel('Forecast range (UTC)',fontsize=19)
# ax.yaxis.set_label_coords(-0.10, .5) # place of labels
# ax.xaxis.set_label_coords(0.5, -0.14)
# ax.yaxis.set_tick_params(labelsize=16)
# ax.xaxis.set_tick_params(labelsize=16)

# ax.plot(SST_GN84.ech,GKCJ,color=color_GKCJ,label='REF.noOML',linewidth=3)
# # ax.plot(SST_GN84.ech,SST_GL55.SST,color=color_GL55,label='L.26',linewidth=3)
# ax.plot(SST_GN84.ech,SST_GN84.SST,color=color_GN84,label='Xp.OML',linewidth=3)
# ax.plot(SST_GN84.ech,SST_GO4A.SST,color=color_GO4A,label='Xp.OML.buyos',linewidth=3)
# ax.plot(SST_GN84.ech,SST_GO4A_an.SST,color=color_GO4A,label='Xp.OML.buoys.Analysis',linewidth=3,linestyle='--')
# # ax.plot(SST_GN84.ech,SST_GOJQ.SST,color=color_GOJQ,label='R10x4.buoys',linewidth=3)
# ax.plot(SST_GN84.ech,df_diurne_oc.SST,color=color_GN3C,label='SST.Mercator',linewidth=3)
# ax.plot(SST_GN84.ech,df_diurne_bouees.SST,color='lightskyblue',label='Obs.buyos', linewidth=3)
# plt.grid(alpha=1)
# ax.legend(bbox_to_anchor=(0.5, 0., 0.5, 0.8))
# plt.legend(fontsize=18, bbox_to_anchor=(1.03, 0.70))

# # fig.suptitle('SST evolution over the forecast period - 102h - of the '+str(Date_debut)[0:10]+'(lon:'+str(lon_expe)+' lat:'+str(lat_expe)+')',fontsize=16,y=0.96)
# plt.show()

# figname = Dir_im+str(Date_debut)[0:10]+'_lonlat_'+str(lon_expe)+'_'+str(lat_expe)+'_SST_evol_forecast_lessplot.png'
# fig.savefig(figname,dpi=250, format='png',bbox_inches='tight')


#%%plot sans bouees
SST_GKCJ = df_diurne.loc[((df_diurne['expe']=='GKCJ') )][['date','expe','var','level','ech', 'SST']]
GKCJ=SST_GKCJ.SST.iloc[1]
GKCJ=np.linspace(GKCJ,GKCJ,17)

SST_GN84 = df_diurne.loc[((df_diurne['expe']=='GN84'))][['date','expe','var','level','ech', 'SST']] 
#SST_GL55 = df_diurne.loc[((df_diurne['expe']=='GL55'))][['date','expe','var','level','ech', 'SST']] 
SST_GO4A = df_diurne.loc[((df_diurne['expe']=='GO4A'))][['date','expe','var','level','ech', 'SST']] 
SST_GOJQ = df_diurne.loc[((df_diurne['expe']=='GOJQ'))][['date','expe','var','level','ech', 'SST']] 

SST_GO4A_an = df_diurne_an.loc[((df_diurne_an['expe']=='GO4A'))][['date','expe','var','level','ech', 'SST']]

fig, ax = plt.subplots(nrows = 1,ncols = 1, sharex = False, figsize=(8, 5),dpi=150) 

## BIAIS
# ax.set_xticks([])
ax.set_ylabel('SST (°C)',fontsize=19)
# xlabels = df_diurne_atm.ech[0:len(df_diurne_atm.ech):N]
# ax.set_xticks(xind,labels=xlabels,fontsize=15)

ax.set_xlabel('Forecast range (UTC)',fontsize=19)
ax.yaxis.set_label_coords(-0.12, .5) # place of labels
ax.xaxis.set_label_coords(0.5, -0.14)
ax.yaxis.set_tick_params(labelsize=16)
ax.xaxis.set_tick_params(labelsize=16)

ax.plot(SST_GN84.ech,GKCJ,color=color_GKCJ,label='REF.noOML',linewidth=3)
# ax.plot(SST_GN84.ech,SST_GL55.SST,color=color_GL55,label='L.26',linewidth=3)
ax.plot(SST_GN84.ech,SST_GN84.SST,color=color_GN84,label='Xp.OML',linewidth=3)
ax.plot(SST_GN84.ech,SST_GO4A.SST,color=color_GO4A,label='Xp.OML.buyos',linewidth=3)
ax.plot(SST_GN84.ech,SST_GO4A_an.SST,color=color_GO4A,label='Xp.OML.buoys.Analysis',linewidth=3,linestyle='--')
# ax.plot(SST_GN84.ech,SST_GOJQ.SST,color=color_GOJQ,label='R10x4.buoys',linewidth=3)
ax.plot(SST_GN84.ech,df_diurne_oc.SST,color=color_GN3C,label='SST.Mercator',linewidth=3)
# ax.plot(SST_GN84.ech,df_diurne_bouees.SST,color='lightskyblue',label='Obs.buyos', linewidth=3)
plt.grid(alpha=1)
ax.legend(bbox_to_anchor=(0.5, 0., 0.5, 0.8))
plt.legend(fontsize=18, bbox_to_anchor=(0.99, 0.70))

# fig.suptitle('SST evolution over the forecast period - 102h - of the '+str(Date_debut)[0:10]+'(lon:'+str(lon_expe)+' lat:'+str(lat_expe)+')',fontsize=16,y=0.96)
plt.show()

# figname = Dir_im+str(Date_debut)[0:10]+'_lonlat_'+str(lon_expe)+'_'+str(lat_expe)+'_SST_evol_forecast_nobuoys.png'
# fig.savefig(figname,dpi=250, format='png',bbox_inches='tight')

#%%

SST_GKCJ = df_diurne.loc[((df_diurne['expe']=='GKCJ') )][['date','expe','var','level','ech', 'SST']]
GKCJ=SST_GKCJ.SST.iloc[1]
GKCJ=np.linspace(GKCJ,GKCJ,17)

SST_GN84 = df_diurne.loc[((df_diurne['expe']=='GN84'))][['date','expe','var','level','ech', 'SST']] 
#SST_GL55 = df_diurne.loc[((df_diurne['expe']=='GL55'))][['date','expe','var','level','ech', 'SST']] 
# SST_GO4A = df_diurne.loc[((df_diurne['expe']=='GO4A'))][['date','expe','var','level','ech', 'SST']] 
# SST_GOJQ = df_diurne.loc[((df_diurne['expe']=='GOJQ'))][['date','expe','var','level','ech', 'SST']] 
SST_GPPS = df_diurne.loc[((df_diurne['expe']=='GPPS'))][['date','expe','var','level','ech', 'SST']] 
SST_GPPU = df_diurne.loc[((df_diurne['expe']=='GPPU'))][['date','expe','var','level','ech', 'SST']] 
SST_GPPN = df_diurne.loc[((df_diurne['expe']=='GPPN'))][['date','expe','var','level','ech', 'SST']] 

SST_GNBA = df_diurne.loc[((df_diurne['expe']=='GNBA'))][['date','expe','var','level','ech', 'SST']] 
# SST_GO4A_an = df_diurne_an.loc[((df_diurne_an['expe']=='GO4A'))][['date','expe','var','level','ech', 'SST']]

fig, ax = plt.subplots(nrows = 1,ncols = 1, sharex = False, figsize=(8, 5),dpi=150) 

## BIAIS
# ax.set_xticks([])
ax.set_ylabel('SST (°C)',fontsize=19)
# xlabels = df_diurne_atm.ech[0:len(df_diurne_atm.ech):N]
# ax.set_xticks(xind,labels=xlabels,fontsize=15)

ax.set_xlabel('Forecast range (UTC)',fontsize=19)
ax.yaxis.set_label_coords(-0.12, .5) # place of labels
ax.xaxis.set_label_coords(0.5, -0.14)
ax.yaxis.set_tick_params(labelsize=16)
ax.xaxis.set_tick_params(labelsize=16)

# ax.plot(SST_GN84.ech,GKCJ,color=color_GKCJ,label='REF.noOML',linewidth=3)

# ax.plot(SST_GN84.ech,SST_GL55.SST,color=color_GL55,label='L.26',linewidth=3)
ax.plot(SST_GN84.ech,SST_GN84.SST,color=color_GN84,label='Xp.OML',linewidth=3)
ax.plot(SST_GN84.ech,SST_GNBA.SST,color='mediumseagreen',label='Xp.turbidite_GNBA',linewidth=3)
# ax.plot(SST_GN84.ech,SST_GPPS.SST,color='blue',label='Xp.TKEx10_GPPS',linewidth=3,linestyle='--')
# ax.plot(SST_GN84.ech,SST_GPPN.SST,color='skyblue',label='Xp.TKEdiv10_GPPS',linewidth=3,linestyle='--')

# ax.plot(SST_GN84.ech,SST_GPPU.SST,color='tomato',label='Xp.coeffdiv5_GPPU',linewidth=3,linestyle='--')


# ax.plot(SST_GN84.ech,SST_GO4A.SST,color=color_GO4A,label='Xp.OML.buyos',linewidth=3)
# ax.plot(SST_GN84.ech,SST_GO4A_an.SST,color=color_GO4A,label='Xp.OML.buoys.Analysis',linewidth=3,linestyle='--')
# ax.plot(SST_GN84.ech,SST_GOJQ.SST,color=color_GOJQ,label='R10x4.buoys',linewidth=3)
# ax.plot(SST_GN84.ech,df_diurne_oc.SST,color=color_GN3C,label='SST.Mercator',linewidth=3)

# ax.plot(SST_GN84.ech,df_diurne_bouees.SST,color='lightskyblue',label='Obs.buyos', linewidth=3)
plt.grid(alpha=1)
ax.legend(bbox_to_anchor=(0.5, 0., 0.5, 0.8))
plt.legend(fontsize=18, bbox_to_anchor=(0.99, 0.70))



figname = Dir_im+str(Date_debut)[0:10]+'_lonlat_'+str(lon_expe)+'_'+str(lat_expe)+'_SST_evol_test_coeff_donnees_GNBA.png'
fig.savefig(figname,dpi=250, format='png',bbox_inches='tight')
# fig.suptitle('SST evolution over the forecast period - 102h - of the '+str(Date_debut)[0:10]+'(lon:'+str(lon_expe)+' lat:'+str(lat_expe)+')',fontsize=16,y=0.96)
plt.show()
#%%
print('lon_bouees:',lon_bouees_pt)
print('lat_bouees:',lat_bouees_pt)
print('sst bouees mean:',np.round(df_diurne_bouees.SST.mean(),3))
print('sst bouees min:',np.round(df_diurne_bouees.SST.min(),3))
print('sst bouees max:',np.round(df_diurne_bouees.SST.max(),3))
print('sst bouees std:',np.round(df_diurne_bouees.SST.std(),3))


print('sst merc mean:',np.round(df_diurne_oc.SST.mean(),3))
print('sst merc min:',np.round(df_diurne_oc.SST.min(),3))
print('sst merc max:',np.round(df_diurne_oc.SST.max(),3))
print('sst merc std:',np.round(df_diurne_oc.SST.std(),3))

print('sst forecast Xp.finale.buyos mean:',np.round(SST_GO4A.SST.mean(),3))
print('sst forecast Xp.finale.buyos min:',np.round(SST_GO4A.SST.min(),3))
print('sst forecast Xp.finale.buyos max:',np.round(SST_GO4A.SST.max(),3))
print('sst forecast Xp.finale.nuyos std:',np.round(SST_GO4A.SST.std(),3))

print('sst analyse Xp.finale.buyos mean:',np.round(SST_GO4A_an.SST.mean(),3))
print('sst analyse Xp.finale.buyos min:',np.round(SST_GO4A_an.SST.min(),3))
print('sst analyse Xp.finale.buyos max:',np.round(SST_GO4A_an.SST.max(),3))
print('sst analyse Xp.finale.buyos std:',np.round(SST_GO4A_an.SST.std(),3))

print('sst double:',np.round(SST_GKCJ.SST.mean(),3))
