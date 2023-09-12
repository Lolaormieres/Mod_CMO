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
#MINMAX=[0,1]
#MINMAX=[270,280] #rajout
#MINMAX=[[0,1],[270,280],[0,4.],[-40,1]]
# MINMAX=[-0.5,0.5,-0.5,0.5]
# MINMAX=[-2.5,2.5]
MINMAX=[-2,2]
# MINMAX=[[-2,2]] # ,[-10,10],[-10,10],[-10,10]]
#MINMAX=[[-10,10]]#,[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]]
#MINMAX=[-2.3,10]

name=['SFX.TEMP_OC1']#,'SFX.TEMP_OC12']
      #'SFX.TEMP_OC2','SFX.TEMP_OC3','SFX.TEMP_OC4','SFX.TEMP_OC5','SFX.TEMP_OC6','SFX.TEMP_OC7','SFX.TEMP_OC8','SFX.TEMP_OC9','SFX.TEMP_OC10','SFX.TEMP_OC11','SFX.TEMP_OC12']#,'SFX.TEMP_OC13','SFX.TEMP_OC14','SFX.TEMP_OC15','SFX.TEMP_OC16','SFX.TEMP_OC17','SFX.TEMP_OC18','SFX.TEMP_OC19','SFX.TEMP_OC20','SFX.TEMP_OC21','SFX.TEMP_OC22','SFX.TEMP_OC23','SFX.TEMP_OC24','SFX.TEMP_OC25','SFX.TEMP_OC26']
name_sal=['SFX.SALT_OC1']#,'SFX.SALT_OC2','SFX.SALT_OC3','SFX.SALT_OC4','SFX.SALT_OC5','SFX.SALT_OC6','SFX.SALT_OC7','SFX.SALT_OC8','SFX.SALT_OC9','SFX.SALT_OC10','SFX.SALT_OC11','SFX.SALT_OC12']#,'SFX.SALT_OC13','SFX.SALT_OC14','SFX.SALT_OC15','SFX.SALT_OC16','SFX.SALT_OC17','SFX.SALT_OC18','SFX.SALT_OC19','SFX.SALT_OC20''SFX.SALT_OC21','SFX.SALT_OC22','SFX.SALT_OC23','SFX.SALT_OC24','SFX.SALT_OC25','SFX.SALT_OC26']
#name_mer=['SURFSST.CLIM.']
name_mer=['SURF.THETAO1']#,'SURF.THETAO12','SURF.THETAO2','SURF.THETAO3','SURF.THETAO4','SURF.THETAO5','SURF.THETAO6','SURF.THETAO7','SURF.THETAO8','SURF.THETAO9','SURF.THETAO10','SURF.THETAO11','SURF.THETAO12','SURF.THETAO13','SURF.THETAO14','SURF.THETAO15','SURF.THETAO16','SURF.THETAO17','SURF.THETAO18','SURF.THETAO19','SURF.THETAO20','SURF.THETAO21','SURF.THETAO22','SURF.THETAO23','SURF.THETAO24','SURF.THETAO25','SURF.THETAO26']
name_mer_sal=['SURF.SALINO1','SURF.SALINO2','SURF.SALINO3','SURF.SALINO4','SURF.SALINO5','SURF.SALINO6','SURF.SALINO7','SURF.SALINO8','SURF.SALINO9','SURF.SALINO10','SURF.SALINO11','SURF.SALINO12']#,'SURF.SALINO13','SURF.SALINO14','SURF.SALINO15','SURF.SALINO16','SURF.SALINO17','SURF.SALINO18','SURF.SALINO19','SURF.SALINO20','SURF.SALINO21','SURF.SALINO22','SURF.SALINO23','SURF.SALINO24','SURF.SALINO25','SURF.SALINO26']
name_u=['SFX.UCUR_OC1','SFX.UCUR_OC2','SFX.UCUR_OC3','SFX.UCUR_OC4','SFX.UCUR_OC5','SFX.UCUR_OC6','SFX.UCUR_OC7','SFX.UCUR_OC8','SFX.UCUR_OC9','SFX.UCUR_OC10','SFX.UCUR_OC11','SFX.UCUR_OC12','SFX.UCUR_OC13','SFX.UCUR_OC14','SFX.UCUR_OC15','SFX.UCUR_OC16','SFX.UCUR_OC17','SFX.UCUR_OC18','SFX.UCUR_OC19','SFX.UCUR_OC20','SFX.UCUR_OC21','SFX.UCUR_OC22','SFX.UCUR_OC23','SFX.UCUR_OC24','SFX.UCUR_OC25','SFX.UCUR_OC26']
name_v=['SFX.VCUR_OC1','SFX.VCUR_OC2','SFX.VCUR_OC3','SFX.VCUR_OC4','SFX.VCUR_OC5','SFX.VCUR_OC6','SFX.VCUR_OC7','SFX.VCUR_OC8','SFX.VCUR_OC9','SFX.VCUR_OC10','SFX.VCUR_OC11','SFX.VCUR_OC12','SFX.VCUR_OC13','SFX.VCUR_OC14','SFX.VCUR_OC15','SFX.VCUR_OC16','SFX.VCUR_OC17','SFX.VCUR_OC18','SFX.VCUR_OC19','SFX.VCUR_OC20','SFX.VCUR_OC21','SFX.VCUR_OC22','SFX.VCUR_OC23','SFX.VCUR_OC24','SFX.VCUR_OC25','SFX.VCUR_OC26']
#MINMAX=[[300,310]] # ,[-10,10],[-10,10]]
#MINMAX=[[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300],[274,300]]

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




zone='Eurat'
expe='GN84'
# expe='GO48'


# expe='GN3C'

xmer ='Mercator'
#xpref='GKOR'
expe_drown='GMPE'
#geometry='global798c22'
geometry='global1798'
#term='24'
#chemin='/d0/images/ormieresl/'
chemin='/home/ormieresl/Scripts/Cartes'
chemin_pdg='/home/ormieresl/'

Dir_im='/cnrm/recyf/Data/users/ormieresl/plot_cartes_diurne/'

## ARP
cfields=''
cblock='surfan'
ckind='analysis'
cut='assim'
#cblock='forecast'
#ckind='historic'
cblock_for='forecast'
ckind_for='historic'
cmodel='surfex'
cterm_0='0'
cmodel_for='surfex'
cut_for='production'
cterm_6='6'
cterm_12='12'
cterm_18='18'
cterm_24='24'
cterm_30='30'
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
cmodel_oc='arpege'
cut_oc='assim'
cterm_oc=''


long=170
lati=-81

#domaines lonmin,lonmax,latmin,matmax
#zones = {"nordat": [-80,0,20,70], "eurat": [-35,45,20,72], "tropiques": [0,360,-20,20],"hn20":[0,360,20,90],"hs20":[0,360,-90,-20],"glob":[0,360,-90,90],"med": [-3,16,30,44]}

def get_file(xpid,date,rep,typfic,fields,mod,chain,ech):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geometry,'format':'fa',
          'local':'tmpcartes_moy.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr','filling':'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpcartes_moy.fic',openmode='r')
    os.remove('tmpcartes_moy.fic')
    return r


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
        
        #1) CMO
        crs=ccrs.PlateCarree()
        fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='biais_v4_grey')
        # ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
        # ax.title.set_text(name[i]+'_'+str(expe)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-')
        # ax.title.set_text(name[i]+'_'+str(expe)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours '+str(cterm_18)+'-'+str(cterm_6)+'-GLOB-')
        # ax.title.set_text('Diurnal SST range in OML model (xp.final), average over september (18UTC - 6UTC)', fontsize=26)
        # clb=fig.colorbar
        # ax.tick_params(axis='y', labelsize=15)
        # clb.ax.set_title('Bias (°C)',fontsize=19)
        plt.yticks(fontsize=20)  
        
        ax.set_title('Diurnal SST range in OML model average over september (18UTC - 6UTC), expe:'+str(expe), fontsize=26,y=1.02)
       
        
        # ax.set_yticks(np.arange(20,71,10), crs=ccrs.PlateCarree())
        # ax.yaxis.set_major_formatter(lat_formatter)
        ax.tick_params(axis='y', labelsize=18)
        #Define the xticks for latitude
        # ax.set_xticks(np.arange(-35,46,10), crs=ccrs.PlateCarree())
        
        ax.tick_params(axis='y', labelsize=18)
        ax.tick_params(axis='x', labelsize=18)
        plt.show()
        
        # fig.savefig(chemin+name[i]+'_'+str(expe)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-'+cblock+'-'+ckind+str(cterm)+'.png')
        # fig.colorbar(fig, spacing='uniform', orientation='vertical')
        figname = Dir_im+str(expe)+'_'+str(zone)+'_'+str(nb)+'jours.png'
        fig.savefig(figname,dpi=300, format='png',bbox_inches='tight')
        figname = Dir_im+str(expe)+'_'+str(zone)+'_'+str(nb)+'jours.pdf'
        fig.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')
      
 
        # # 2) MERCATOR
        # crs=ccrs.PlateCarree()
        # fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='biais_v4_grey')
        # # ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        # ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
        # # ax.title.set_text(name[i]+'_'+str(expe)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-')
        # # ax.title.set_text(name[i]+'_'+str(expe)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours '+str(cterm_18)+'-'+str(cterm_6)+'-GLOB-')
        # # ax.title.set_text('Diurnal SST range in OML model (xp.final), average over september (18UTC - 6UTC)', fontsize=26)
        # # clb=fig.colorbar
        # # ax.tick_params(axis='y', labelsize=15)
        # # clb.ax.set_title('Bias (°C)',fontsize=19)
        # plt.yticks(fontsize=20)  
        
        # ax.set_title('Diurnal SST range in the REF. instantaneous SST Mercator, average over september (18UTC - 6UTC)', fontsize=26,y=1.02)
       
        
        # # ax.set_yticks(np.arange(20,71,10), crs=ccrs.PlateCarree())
        # # ax.yaxis.set_major_formatter(lat_formatter)
        # ax.tick_params(axis='y', labelsize=18)
        # #Define the xticks for latitude
        # # ax.set_xticks(np.arange(-35,46,10), crs=ccrs.PlateCarree())
        
        # ax.tick_params(axis='y', labelsize=18)
        # ax.tick_params(axis='x', labelsize=18)
        # plt.show()
        
        # # fig.savefig(chemin+name[i]+'_'+str(expe)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-'+cblock+'-'+ckind+str(cterm)+'.png')
        # # fig.colorbar(fig, spacing='uniform', orientation='vertical')
        # figname = Dir_im+str(expe)+'_'+str(zone)+'_diurnal_amplitude_mercator_'+str(nb)+'jours.png'
        # fig.savefig(figname,dpi=300, format='png',bbox_inches='tight')
        # figname = Dir_im+str(expe)+'_'+str(zone)+'_diurnal_amplitude_mercatir_'+str(nb)+'jours.pdf'
        # fig.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')
        
        
        
        
        
        # crs=ccrs.PlateCarree()
        # fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='biais_v4_grey')
        # # ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        # ax.set_extent([-35, 45, 20, 72], ccrs.PlateCarree())
        # # ax.title.set_text(name[i]+'_'+str(expe)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-')
        # # ax.title.set_text(name[i]+'_'+str(expe)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours '+str(cterm_18)+'-'+str(cterm_6)+'-GLOB-')
        
        # ax.tick_params(axis='y', labelsize=18)
        # ax.title.set_text('Diurnal SST range in the REF. instantaneous SST Mercator, average over september (18UTC - 6UTC)')
        
        # ax.set_title('Diurnal SST range in the REF. instantaneous SST Mercator, average over september (18UTC - 6UTC)', fontsize=26,y=1.02)
        # plt.yticks(fontsize=20)
        # plt.show()
        
        # # fig.savefig(chemin+name[i]+'_'+str(expe)+'-'+str(xmer)+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'jours'+'-GLOB-'+cblock+'-'+ckind+str(cterm)+'.png')
        # # fig.colorbar(fig, spacing='uniform', orientation='vertical')
        # figname = '/cnrm/recyf/Data/users/ormieresl/mean_diurnal_amplitude_mercator'+'_'+str(zone)+'.png'
        # fig.savefig(figname,dpi=300, format='png',bbox_inches='tight')
        # figname = '/cnrm/recyf/Data/users/ormieresl/mean_diurnal_amplitude_mercator'+'_'+str(zone)+'.pdf'
        # fig.savefig(figname,dpi=300, format='pdf',bbox_inches='tight')
       
## Diff toutes les 6 heures avec GN3C
# nb=60
# date_init=datetime(2022,9,2,0)
# date=date_init
# dt = timedelta(hours=6)

# expe=['GL55', 'GKOR', 'GKPH','GM6E', 'GMOT']
# for e in range(len(expe)): 

    



# Eurat
# Long=-35;45,Lati=20;72
# Long=-35, Lati=20 

# # Diff deux exp
# for t in range(nb):
#     print(date)
    
#     try:
#         r = get_file(expe, date, cblock, ckind, cfields, cmodel, cut, cterm)
#         # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
#         # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
#         r_ref = get_file(xpref, date, cblock, ckind, cfields, cmodel, cut, cterm)
#         if t == 0:
#             x_diff = r.readfield(name[0])-(r_ref.readfield(name[0]))
#         else:
#             x_diff=r.readfield(name[0])-(r_ref.readfield(name[0]))+x_diff
#     except:
#         print('Pas de fichier')
#     date=date+dt
    
# x_diff.setdata(x_diff.getdata()/nb)
# plot_diffmean()


# =============================================================================
# CYCLE DIRUNE
# =============================================================================
#-------------------------------
#calcul diff
#-------------------------------
#Diff journaliere
nb=30
date_init=datetime(2022,9,1,0) # put 0 for CMO and 6 for Mercator
date=date_init

dh=timedelta(hours=12)
date_plus=date+dh
dt=timedelta(days=1)

# 1) XP CMO
# TEST : cterm_18, cterm_30
# cterm_0, cterm_12
 
for t in range(nb):
    print(date)
    
#try:
    # r = get_file(expe, date, cblock, ckind, cfields, cmodel, cut, cterm)
    
    # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
    r_moins=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_18)
    r_plus=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_30)
    # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
    # r_mer_moins=get_file(expe,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc)
    # r_mer_plus=get_file(expe,date_plus,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm_oc)
    
    X=r_moins.geometry.get_lonlat_grid()  #lon et lat
    lon=X[0]
    lat=X[1]  
      
    mask=(lon<45)&(lon>-35)&(lat>20)&(lat<72)
      
    if t == 0:
        x_diff = r_moins.readfield(name[0])-(r_plus.readfield(name[0])) 
    else:
        x_diff=r_moins.readfield(name[0])-(r_plus.readfield(name[0]))+x_diff
    # except:
    #     print('Pas de fichier')
    date=date+dt
    date_plus=date+dh


print('mean diurnal amplitudes total - CMO: ',np.round(x_diff.data.mean(),2))
print('min diurnal amplitudes total - CMO: ',np.round(x_diff.data.min(),2))
print('max diurnal amplitudes total - CMO: ',np.round(x_diff.data.max(),2))

print('mean diurnal amplitudes - CMO: ',np.round(x_diff.data[mask].mean(),2))
print('min diurnal amplitudes - CMO: ',np.round(x_diff.data[mask].min(),2))
print('max diurnal amplitudes - CMO: ',np.round(x_diff.data[mask].max(),2))

x_diff.setdata(x_diff.getdata()/nb)

print('mean diurnal amplitudes total - CMO: ',np.round(x_diff.data.mean(),2))
print('min diurnal amplitudes total - CMO: ',np.round(x_diff.data.min(),2))
print('max diurnal amplitudes total - CMO: ',np.round(x_diff.data.max(),2))

print('mean diurnal amplitudes - CMO: ',np.round(x_diff.data[mask].mean(),2))
print('min diurnal amplitudes - CMO: ',np.round(x_diff.data[mask].min(),2))
print('max diurnal amplitudes - CMO: ',np.round(x_diff.data[mask].max(),2))
plot_diffmean()



# ## 2) XP MERCATOR
# TEST : cterm_18, cterm_30
# cterm_0, cterm_12

# expe=['GL55','GKOR','GKPH','GN84','GO4A','GNYG']
# for e in range(len(expe)):
#     print(e, expe[e])
# Diff avec MERCATOR    
# for t in range(nb):
#     print(date)

# # try:
#     # r = get_file(expe, date, cblock, ckind, cfields, cmodel, cut, cterm)

#     # r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
#     # r_moins=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_6)
#     # r_plus=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_18)
#     # r_mer = get_file(expe, date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
#     r_mer_moins = get_file(expe, date, cblock_oc, ckind_oc,
#                             cfields_oc, cmodel_oc, cut_oc, cterm_oc)
#     r_mer_plus = get_file(expe, date_plus, cblock_oc,
#                           ckind_oc, cfields_oc, cmodel_oc, cut_oc, cterm_oc)

#     Y = r_mer_moins.geometry.get_lonlat_grid()  # lon et lat
#     lon_mer = Y[0]
#     lat_mer = Y[1]
#     f = r_mer_moins.readfield(name_mer[0])

#     mask = (lon_mer < 45) & (lon_mer > -35) & (lat_mer >
#                                                 20) & (lat_mer < 72) & (f.data > 200) & (f.data < 400)

#     if t == 0:
#         x_diff = r_mer_plus.readfield(
#             name_mer[0])-(r_mer_moins.readfield(name_mer[0]))

#     else:
#         x_diff = r_mer_plus.readfield(
#             name_mer[0])-(r_mer_moins.readfield(name_mer[0]))+x_diff

#     # except:
#     #     print('Pas de fichier')
#     date = date+dt
#     date_plus = date+dh


# print('mean diurnal amplitudes total - MERC: ', np.round(x_diff.data.mean(), 2))
# print('min diurnal amplitudes total - MERC: ', np.round(x_diff.data.min(), 2))
# print('max diurnal amplitudes total - MERC: ', np.round(x_diff.data.max(), 2))

# print('mean diurnal amplitudes - MERC: ',
#       np.round(x_diff.data[mask].mean(), 2))
# print('min diurnal amplitudes - MERC: ', np.round(x_diff.data[mask].min(), 2))
# print('max diurnal amplitudes - MERC: ', np.round(x_diff.data[mask].max(), 2))


# x_diff.setdata(x_diff.getdata()/nb)


# print('mean diurnal amplitudes total - MERC: ', np.round(x_diff.data.mean(), 2))
# print('min diurnal amplitudes total - MERC: ', np.round(x_diff.data.min(), 2))
# print('max diurnal amplitudes total - MERC: ', np.round(x_diff.data.max(), 2))

# print('mean diurnal amplitudes - MERC: ',
#       np.round(x_diff.data[mask].mean(), 2))
# print('min diurnal amplitudes - MERC: ', np.round(x_diff.data[mask].min(), 2))
# print('max diurnal amplitudes - MERC: ', np.round(x_diff.data[mask].max(), 2))

# plot_diffmean()
