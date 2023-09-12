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
MINMAX=[-5,5]
#MINMAX=[[-10,10]]#,[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]]
#MINMAX=[-2.3,10]

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

# Lecture des niveaux de la CMO
filename = 'level2.txt'
level = np.loadtxt(filename, delimiter=',', dtype=float)
level_12=[]
level_12.append(level[0:12])
# level_19=[]
# level_19.append(level[0:19])
# level_26=[]
# level_26.append(level[0:26])
# level_21=[]
# level_21.append(level[0:21])
level_12=np.array(level_12)
level_12=level_12.T
# level_26=np.array(level_26)
# level_26=level_26.T
# level_21=np.array(level_21)
# level_21=level_21.T

level=level_12

expe='GL55'
expe2='GL55'
expe3='GKPH'
xpref='MER'

xmer ='Mercator'
#xpref='GKOR'
expe_drown='GMPE'
#geometry='global798c22'
geometry='global1798'
#term='24'
#chemin='/d0/images/ormieresl/'
chemin='/home/ormieresl/Scripts/'
chemin_pdg='/home/ormieresl/'

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
cterm='0'
cmodel_for='surfex'
cut_for='production'
cterm_6='6'
cterm_12='12'
cterm_18='18'
cterm_24='24'
cterm_48='48'
cterm_54='54'
cterm_60='60'
cterm_p6='6'
cterm_for='96'
## Mer
#cblock_oc='c933_iter1_10_testITM'
cblock_oc='c933_iter1_10_testinverse'
cblock_oc='c933'
ckind_oc='geofields'
cfields_oc='ocean'
#cfields_oc='sst'
cmodel_oc='arpege'
cut_oc='assim'
long=170
lati=-81

dir_fig = '/cnrm/recyf/Data/users/ormieresl/plot_sst_cyclones/'
def get_file(xpid,date,rep,typfic,fields,mod,chain,ech):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geometry,'format':'fa',
          'local':'tmpcartes.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr','filling':'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmpcartes.fic',openmode='r')
    os.remove('tmpcartes.fic')
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


def plot_diffmean():
    for i in range(len(name)):
        print(i,name[i])
        crs=ccrs.PlateCarree()
        fig, ax = x_diff.cartoplot(projection=crs,minmax=MINMAX,plot_method='scatter',colormap='biais_v4_grey_v3')
        # ax.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        ax.set_extent([-100, -70, 15, 35], ccrs.PlateCarree())

        # # taken fromd undef_zoom tickbar --> make bigger lon/lat
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = None
        gl.ylabels_left = True
        gl.xlabels_bottom = True
        gl.ylabels_right = None
        gl.xlabel_style = {'size': 17, 'color': 'black', 'weight':'bold'}
        gl.ylabel_style = {'size': 17, 'color': 'black','weight':'bold'}
        
        ax.tick_params(axis='y', labelsize=18)
        plt.yticks(fontsize=20) 
        ax.tick_params(axis='x', labelsize=18)
        
        
        #ax.title.set_text(name[i]+'_'+str(expe[e])+
                         # '_'+str(date)[0:10]+'-'+str(date)[11:13])
        
        plt.suptitle(str(name[0])+'_'+str(expe)+'_Date: '+str(date)
             [0:10]+'-'+str(date)[11:13]+'_P96-P0', fontsize=26)
        plt.show()
        fig.savefig(dir_fig+name[i]+'_'+str(expe[e])+'_'+str(date_init)[0:10]+'-'+str(date_init)[11:13]+'_'+str(nb)+'-IAN-'+'P'+str(cterm_for)+'.png')
# r.close()
# r_ref.close()
#-------------------------------
#calcul diff
#-------------------------------
#Diff journaliere
nb=1
date_init=datetime(2022,9,28,0)
date=date_init
dt=timedelta(days=1)


## Diff toutes les 6 heures avec GN3C
# nb=60
# date_init=datetime(2022,9,2,0)
# date=date_init
# dt = timedelta(hours=6)

# 


# for t in range(nb):
#     print(date)
#     ## Enleve GMOT 
#     expe=['GN84']
#     for e in range(len(expe)): 
        
        
#         # r = get_file(expe[e], date, cblock, ckind, cfields, cmodel, cut, cterm)
#         r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
#         r_mer = get_file(expe[e], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
#         if t == 0:
#             x_diff = r_for.readfield(name[0])-(r_mer.readfield(name_mer[0])-273.15)
#         else:
#             x_diff=r_for.readfield(name[0])-(r_mer.readfield(name_mer[0])-273.15)+x_diff
        
#         # date=date+dt
    
# # x_diff.setdata(x_diff.getdata()/nb)
#         plot_diffmean() 

expe_ref='GP7E'
for t in range(nb):
    print(date)
    ## Enleve GMOT 
    expe=['GPGP']
    for e in range(len(expe)): 
        
        
        # r = get_file(expe[e], date, cblock, ckind, cfields, cmodel, cut, cterm)
        r_for_moins=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm)
        r_for_plus=get_file(expe[e],date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
        
        r_for_moins_ref=get_file(expe_ref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm)
        r_for_plus_ref=get_file(expe_ref,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
        # r_mer = get_file(expe[e], date, cblock_oc, ckind_oc,cfields_oc, cmodel_oc, cut_oc, cterm)
       
        x_diff = (r_for_plus.readfield(name[0])-(r_for_moins.readfield(name[0]))) 
        #- (r_for_plus_ref.readfield(name[0])-(r_for_moins_ref.readfield(name[0]))) 
         
        x_diff = (r_for_plus.readfield(name[0])-(r_for_plus_ref.readfield(name[0]))) 
        
        # date=date+dt
    
# x_diff.setdata(x_diff.getdata()/nb)
        plot_diffmean() 