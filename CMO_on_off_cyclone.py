#! /usr/bin/python

# coding: utf-8

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
epygram.init_env()




#name=['SFX.SIC','SFX.SST','SFX.ICEHSI_1','SFX.ICETSF_1']
#name=['SFX.SIC']
#name=['SFX.SST']

# MINMAX=[-2.5,2.5]
MINMAX=[-5,5]
MINMAX_values=[0,100]


name=['SFX.SST']
# name=['SFX.TEMP_OC1']#,'SFX.TEMP_OC12']
      #'SFX.TEMP_OC2','SFX.TEMP_OC3','SFX.TEMP_OC4','SFX.TEMP_OC5','SFX.TEMP_OC6','SFX.TEMP_OC7','SFX.TEMP_OC8','SFX.TEMP_OC9','SFX.TEMP_OC10','SFX.TEMP_OC11','SFX.TEMP_OC12']#,'SFX.TEMP_OC13','SFX.TEMP_OC14','SFX.TEMP_OC15','SFX.TEMP_OC16','SFX.TEMP_OC17','SFX.TEMP_OC18','SFX.TEMP_OC19','SFX.TEMP_OC20','SFX.TEMP_OC21','SFX.TEMP_OC22','SFX.TEMP_OC23','SFX.TEMP_OC24','SFX.TEMP_OC25','SFX.TEMP_OC26']
name_sal=['SFX.SALT_OC1']#,'SFX.SALT_OC2','SFX.SALT_OC3','SFX.SALT_OC4','SFX.SALT_OC5','SFX.SALT_OC6','SFX.SALT_OC7','SFX.SALT_OC8','SFX.SALT_OC9','SFX.SALT_OC10','SFX.SALT_OC11','SFX.SALT_OC12']#,'SFX.SALT_OC13','SFX.SALT_OC14','SFX.SALT_OC15','SFX.SALT_OC16','SFX.SALT_OC17','SFX.SALT_OC18','SFX.SALT_OC19','SFX.SALT_OC20''SFX.SALT_OC21','SFX.SALT_OC22','SFX.SALT_OC23','SFX.SALT_OC24','SFX.SALT_OC25','SFX.SALT_OC26']
#name_mer=['SURFSST.CLIM.']
name_mer=['SURF.THETAO1']#,'SURF.THETAO12','SURF.THETAO2','SURF.THETAO3','SURF.THETAO4','SURF.THETAO5','SURF.THETAO6','SURF.THETAO7','SURF.THETAO8','SURF.THETAO9','SURF.THETAO10','SURF.THETAO11','SURF.THETAO12','SURF.THETAO13','SURF.THETAO14','SURF.THETAO15','SURF.THETAO16','SURF.THETAO17','SURF.THETAO18','SURF.THETAO19','SURF.THETAO20','SURF.THETAO21','SURF.THETAO22','SURF.THETAO23','SURF.THETAO24','SURF.THETAO25','SURF.THETAO26']
name_mer_sal=['SURF.SALINO1','SURF.SALINO2','SURF.SALINO3','SURF.SALINO4','SURF.SALINO5','SURF.SALINO6','SURF.SALINO7','SURF.SALINO8','SURF.SALINO9','SURF.SALINO10','SURF.SALINO11','SURF.SALINO12']#,'SURF.SALINO13','SURF.SALINO14','SURF.SALINO15','SURF.SALINO16','SURF.SALINO17','SURF.SALINO18','SURF.SALINO19','SURF.SALINO20','SURF.SALINO21','SURF.SALINO22','SURF.SALINO23','SURF.SALINO24','SURF.SALINO25','SURF.SALINO26']
name_u=['SFX.UCUR_OC1','SFX.UCUR_OC2','SFX.UCUR_OC3','SFX.UCUR_OC4','SFX.UCUR_OC5','SFX.UCUR_OC6','SFX.UCUR_OC7','SFX.UCUR_OC8','SFX.UCUR_OC9','SFX.UCUR_OC10','SFX.UCUR_OC11','SFX.UCUR_OC12','SFX.UCUR_OC13','SFX.UCUR_OC14','SFX.UCUR_OC15','SFX.UCUR_OC16','SFX.UCUR_OC17','SFX.UCUR_OC18','SFX.UCUR_OC19','SFX.UCUR_OC20','SFX.UCUR_OC21','SFX.UCUR_OC22','SFX.UCUR_OC23','SFX.UCUR_OC24','SFX.UCUR_OC25','SFX.UCUR_OC26']
name_v=['SFX.VCUR_OC1','SFX.VCUR_OC2','SFX.VCUR_OC3','SFX.VCUR_OC4','SFX.VCUR_OC5','SFX.VCUR_OC6','SFX.VCUR_OC7','SFX.VCUR_OC8','SFX.VCUR_OC9','SFX.VCUR_OC10','SFX.VCUR_OC11','SFX.VCUR_OC12','SFX.VCUR_OC13','SFX.VCUR_OC14','SFX.VCUR_OC15','SFX.VCUR_OC16','SFX.VCUR_OC17','SFX.VCUR_OC18','SFX.VCUR_OC19','SFX.VCUR_OC20','SFX.VCUR_OC21','SFX.VCUR_OC22','SFX.VCUR_OC23','SFX.VCUR_OC24','SFX.VCUR_OC25','SFX.VCUR_OC26']
#MINMAX=[[300,310]] # ,[-10,10],[-10,10]]
#MINMAX=[[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300]]

# ATMO
name_atm=[{'shortName':'ssrd'},{'shortName':'10u'},{'shortName':'10v'}] #,{'shortNameECMF':'100v'}] 
# name_atm=[{'shortName':'10u'}] 
# name_atm=[{'shortName':'10v'}] 
#t temp level 0 level 0
# name_atm=[{'shortName':'t'}] 
# name: 'Relative humidity',
# shortName: 'r',
# name_atm=[{'shortName':'r'}] 
# name: '2 metre temperature',
# shortName: '2t',
name_atm=[{'shortNameECMF':'2t'}] 
# name: '2 metre relative humidity',
# shortName: '2r',
# name_atm=[{'shortName':'2r'}] 
name_atm=[{'shortName':'lcc'}] 
short_name=['lcc'] 

# Lecture des niveaux de la CMO
filename = 'level2.txt'
level = np.loadtxt(filename, delimiter=',', dtype=float)
print('level',level)
print(type(level))
level_12=[]
level_12.append(level[0:12])
#print(level_12)
level_19=[]
level_19.append(level[0:19])
level_26=[]
level_26.append(level[0:26])
level_21=[]
level_21.append(level[0:21])



level_12=np.array(level_12)
level_12=level_12.T
level_26=np.array(level_26)
level_26=level_26.T
level_21=np.array(level_21)
level_21=level_21.T
#print('dim',level_26.ndim)
#print('level_12',level_12)
#print('level_26',level_26)

level=level_12

#Diff journaliere
nb=6
date_init=datetime(2022,9,14,0)
date=date_init
dt_mer=timedelta(days=4)
dt=timedelta(days=1)

date_valid=date+dt_mer


expe='GO50'
expe2='GM6E'
expe3='GKPH'
xpref='GN84'

xmer ='Mercator'
#xpref='GKOR'
expe_drown='GMPE'
#geometry='global798c22'
geometry='global1798'
#term='24'
#chemin='/d0/images/ormieresl/'
chemin='/home/ormieresl/Scripts/Cartes'
chemin_pdg='/home/ormieresl/'
dir_images='/home/ormieresl/Routines/Storm_septembre/'

long=170
lati=-81



    ## ARP
cfields=''
cblock='surfan'
ckind='analysis'
cut='assim'
formats='fa'
fill_cmo='surf'
cterm='0'
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
cterm_for='96'
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
# grid.arpege-forecast.glob01+0090:00.grib # Fichier atm qu'on cherche à lire
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

#%%
def get_file(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geo,'format':formats,'filling':fill,
          'local':'tmpmer.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpmer.fic',openmode='r')
    os.remove('tmpmer.fic')
    return r

def get_file_atm(xpid,date,rep,typfic,fields,mod,chain,ech,geo,formats,fill):
    #MODEL
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

### Ceate Dataframe
raw_data = {'date': [''],
                'expe': [''],
                'var': [''],
                'level': [''],
                'biais': [''],
                'lon': [''],
                'lat': ['']}
df = pd.DataFrame(raw_data, columns = ['date','expe','var','level','biais','lon','lat'])


## DIFF EXPE
def plot_diffmean():
    for i in range(len(name)):
        print(i,name[i])
        crs=ccrs.PlateCarree()
        fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='seismic') #'biais_v4_grey'
        # ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
        ax.title.set_text(name[i]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        # ax.title.set_text(name[i]+'_'+str(xpref)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        
        #fig.savefig(chemin+name[i]+'_'+str(expe)+'-'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-'+cblock+'-'+ckind+str(cterm)+'.png')


def plot_values():
    for i in range(len(name)):
        crs=ccrs.PlateCarree()
        fig, ax = x.cartoplot(projection=crs,minmax=MINMAX_values,plot_method='scatter',colormap='Blues') #'biais_v4_grey'(len(name)):
        ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
        # ax.title.set_text(name_atm[0]+'_'+str(xpref)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        ax.title.set_text(short_name[0]+'_'+str(expe)+' '+str(date)[0:10]+'-'+str(date)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')

        # ax.title.set_text(name[i]+'_'+str(xpref)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+str(cterm_for)+'_'+str(nb)+'jours'+'-GLOB-')
        
        fig.savefig(dir_images+short_name[0]+'_'+str(expe)+'-'+'_'+str(date)[0:10]+'-'+str(date)[11:13]+'_'+str(cterm_for)+'.png')
        

# =============================================================================
# #DEUX EXPE
# =============================================================================

# for t in range(nb):
#     print(date)
    
#     try:

#         # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
#         r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
#         r_for_ref=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
          
#         # r_mer = get_file(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
#         if t == 0:
#           x_diff = r_for.readfield(name[0])-(r_for_ref.readfield(name[0]))
          
#         else:
#           x_diff=r_for.readfield(name[0])-(r_for_ref.readfield(name[0]))+x_diff
#     except:
#         print('Pas de fichier')
#     date=date+dt
    
# x_diff.setdata(x_diff.getdata()/nb)
# plot_diffmean()


# =============================================================================
# #MERCATOR  
# =============================================================================
# for t in range(nb):
#     print(date)
    
#     try:
#        # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
        # r_for=get_file_atm(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
        # r_mer = get_file_atm(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
#         if t == 0:
#             x_diff = r_for.readfield(name[0])-(r_mer.readfield(name_mer[0]))
#         else:
#             x_diff=r_for.readfield(name[0])-(r_mer.readfield(name_mer[0]))+x_diff
#     except:
#         print('Pas de fichier')
#     date=date+dt
    
# x_diff.setdata(x_diff.getdata()/nb)
# plot_diffmean()



# =============================================================================
# # ATMOSPHERE
# =============================================================================

# # Diff deux exp
# for t in range(nb):
#     print(date)
    
#     try:
#         # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
        
#         # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
#         # r_for_ref=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
        
#         #ATM:
#         r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#         r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
#           # r_mer = get_file(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
#         # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
#         # r_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel, cut, cterm)
#         if t == 0:
#             x_diff = r_for_atm.readfield(name_atm[0])-(r_for_atm_ref.readfield(name_atm[0]))
#         else:
#             x_diff=r_for_atm.readfield(name_atm[0])-(r_for_atm_ref.readfield(name_atm[0]))+x_diff
#     except:
#         print('Pas de fichier')
#     date=date+dt
    
# x_diff.setdata(x_diff.getdata()/nb)
# plot_diffmean()

 
# =============================================================================
# Values
# =============================================================================
for t in range(nb):
    print(date)
    
    try:
        # r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm,geometry,formats,fill_cmo)
        
        # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
        # r_for_ref=get_file(xpref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for,geometry,formats,fill_cmo)
        
        #ATM:
        r_for_atm=get_file_atm(expe,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
        # r_for_atm_ref=get_file_atm(xpref,date,cblock_atm,ckind_atm,cfields_atm,cmodel_atm,cut_atm,cterm_for,geometry_atm,formats_atm,fill_atm)
          # r_mer = get_file(xpref,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc,geometry,formats,fill_oc)
        # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
        # r_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel, cut, cterm)
       
        x = r_for_atm.readfield(name_atm[0])
        # x = r_for_atm_ref.readfield(name_atm[0])
    except:
        print('Pas de fichier')
    date=date+dt
    
    # x.setdata(x)
    plot_values()
