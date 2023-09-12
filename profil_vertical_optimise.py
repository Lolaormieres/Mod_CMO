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
epygram.init_env()


#name=['SFX.TEMP_OC1','SFX.TEMP_OC2','SFX.TEMP_OC3','SFX.TEMP_OC4','SFX.TEMP_OC5','SFX.TEMP_OC6','SFX.TEMP_OC7','SFX.TEMP_OC8','SFX.TEMP_OC9','SFX.TEMP_OC10','SFX.TEMP_OC11','SFX.TEMP_OC12']
name=['SFX.TEMP_OC1','SFX.TEMP_OC2','SFX.TEMP_OC3','SFX.TEMP_OC4','SFX.TEMP_OC5','SFX.TEMP_OC6','SFX.TEMP_OC7','SFX.TEMP_OC8','SFX.TEMP_OC9','SFX.TEMP_OC10','SFX.TEMP_OC11','SFX.TEMP_OC12']#,'SFX.TEMP_OC13','SFX.TEMP_OC14','SFX.TEMP_OC15','SFX.TEMP_OC16','SFX.TEMP_OC17','SFX.TEMP_OC18','SFX.TEMP_OC19','SFX.TEMP_OC20','SFX.TEMP_OC21','SFX.TEMP_OC22','SFX.TEMP_OC23','SFX.TEMP_OC24','SFX.TEMP_OC25','SFX.TEMP_OC26']
name_sal=['SFX.SALT_OC1','SFX.SALT_OC2','SFX.SALT_OC3','SFX.SALT_OC4','SFX.SALT_OC5','SFX.SALT_OC6','SFX.SALT_OC7','SFX.SALT_OC8','SFX.SALT_OC9','SFX.SALT_OC10','SFX.SALT_OC11','SFX.SALT_OC12','SFX.SALT_OC13','SFX.SALT_OC14','SFX.SALT_OC15','SFX.SALT_OC16','SFX.SALT_OC17','SFX.SALT_OC18','SFX.SALT_OC19','SFX.SALT_OC20''SFX.SALT_OC21','SFX.SALT_OC22','SFX.SALT_OC23','SFX.SALT_OC24','SFX.SALT_OC25','SFX.SALT_OC26']
name_mer=['SURF.THETAO1','SURF.THETAO2','SURF.THETAO3','SURF.THETAO4','SURF.THETAO5','SURF.THETAO6','SURF.THETAO7','SURF.THETAO8','SURF.THETAO9','SURF.THETAO10','SURF.THETAO11','SURF.THETAO12']#,'SURF.THETAO13','SURF.THETAO14','SURF.THETAO15','SURF.THETAO16','SURF.THETAO17','SURF.THETAO18','SURF.THETAO19','SURF.THETAO20','SURF.THETAO21','SURF.THETAO22','SURF.THETAO23','SURF.THETAO24','SURF.THETAO25','SURF.THETAO26']
name_mer_sal=['SURF.SALINO1','SURF.SALINO2','SURF.SALINO3','SURF.SALINO4','SURF.SALINO5','SURF.SALINO6','SURF.SALINO7','SURF.SALINO8','SURF.SALINO9','SURF.SALINO10','SURF.SALINO11','SURF.SALINO12']#,'SURF.SALINO13','SURF.SALINO14','SURF.SALINO15','SURF.SALINO16','SURF.SALINO17','SURF.SALINO18','SURF.SALINO19','SURF.SALINO20','SURF.SALINO21','SURF.SALINO22','SURF.SALINO23','SURF.SALINO24','SURF.SALINO25','SURF.SALINO26']
name_u=['SFX.UCUR_OC1','SFX.UCUR_OC2','SFX.UCUR_OC3','SFX.UCUR_OC4','SFX.UCUR_OC5','SFX.UCUR_OC6','SFX.UCUR_OC7','SFX.UCUR_OC8','SFX.UCUR_OC9','SFX.UCUR_OC10','SFX.UCUR_OC11','SFX.UCUR_OC12']#,'SFX.UCUR_OC13','SFX.UCUR_OC14','SFX.UCUR_OC15','SFX.UCUR_OC16','SFX.UCUR_OC17','SFX.UCUR_OC18','SFX.UCUR_OC19','SFX.UCUR_OC20','SFX.UCUR_OC21','SFX.UCUR_OC22','SFX.UCUR_OC23','SFX.UCUR_OC24','SFX.UCUR_OC25','SFX.UCUR_OC26']
name_v=['SFX.VCUR_OC1','SFX.VCUR_OC2','SFX.VCUR_OC3','SFX.VCUR_OC4','SFX.VCUR_OC5','SFX.VCUR_OC6','SFX.VCUR_OC7','SFX.VCUR_OC8','SFX.VCUR_OC9','SFX.VCUR_OC10','SFX.VCUR_OC11','SFX.VCUR_OC12']#,'SFX.VCUR_OC13','SFX.VCUR_OC14','SFX.VCUR_OC15','SFX.VCUR_OC16','SFX.VCUR_OC17','SFX.VCUR_OC18','SFX.VCUR_OC19','SFX.VCUR_OC20','SFX.VCUR_OC21','SFX.VCUR_OC22','SFX.VCUR_OC23','SFX.VCUR_OC24','SFX.VCUR_OC25','SFX.VCUR_OC26']
name_tke=['SFX.TKE_OC1','SFX.TKE_OC2','SFX.TKE_OC3','SFX.TKE_OC4','SFX.TKE_OC5','SFX.TKE_OC6','SFX.TKE_OC7','SFX.TKE_OC8','SFX.TKE_OC9','SFX.TKE_OC10','SFX.TKE_OC11','SFX.TKE_OC12']

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

expe='GMOT'
expe2='GL55'
expe3='GKPH'
xpref='GKOR'

# expe=['GL55','GKOR','GKPH','GN84','GO4A']
#geometry='global798c22'
geometry='global1798'
#term='24'
#chemin='/d0/images/ormieresl/'
chemin='/home/ormieresl/Scripts/'
Dir_im='/cnrm/recyf/Data/users/ormieresl/plot_prof_vertical/'


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
cterm_for='48'
## Mer
cblock_oc='c933'
ckind_oc='geofields'
cfields_oc='ocean'
cmodel_oc='arpege'
cut_oc='assim'
long=5.13
lati=59.39

# long=6
# lati=59

long=-2
lati=46.7

## Definition color Pale
color_GL55='steelblue' ## GL55 '#648FFF'
color_GKOR='hotpink'  ## GKOR
color_GKPH= 'mediumpurple' ## GKPH '#DC267F'
color_GM6E='paleturquoise'
color_GMOT='peachpuff'
color_GN84='r'
color_GO4A='rosybrown'
color_GNSR='palegreen'
color_GNYG='yellow'

color_expe=color_GMOT
color_mer='gold'

def get_file(xpid,date,rep,typfic,fields,mod,chain,ech):
    #MODEL
    fp = {'experiment':xpid,'kind':typfic,'block':rep,'fields':fields,
          'date':date,'term':ech,'geometry':geometry,'format':'fa',
          'local':'tmp_profil.fic','cutoff':chain,'vapp':'arpege',
          'vconf':'4dvarfr','model':mod,'namespace':'vortex.multi.fr','filling':'surf'}
    fcref = toolbox.input(**fp)[0]
    fcref.get()
    r = epygram.formats.resource(filename='tmp_profil.fic',openmode='r')
    return r

# for e in range(len(expe)):
#     print(e,expe[e])

nb=1
date_init=datetime(2022,9,9,0)
date=date_init
dt=timedelta(days=2)
dt_moins=timedelta(hours=6)
date_valid=date+dt
date_moins=date-dt_moins
cterm_moins='6'


T=[[0,0,0,0,0,0,0]]#Tableau 
T=np.array(T)
x=[]
f_mean=[]
#x=np.array(x)
x_mer=[]
f_mer_mean=[]
x_mer=np.array(x_mer)

x_for=[]
f_for_mean=[]
x_for=np.array(x_for)

    
x_for_moins=[]
f_for_moins_mean=[]
x_for_moins=np.array(x_for_moins)


x_for_init=[]
f_for_init_mean=[]
x_for_init=np.array(x_for_init)

x_mer_init=[]
f_mer_init_mean=[]
x_mer_init=np.array(x_mer_init)

x_anaprod=[]
f_anaprod_mean=[]
x_anaprod=np.array(x_anaprod)  

x_ref=[]
f_ref_mean=[]
x_ref=np.array(x_ref)   

#r_ref=get_file(xpref,date)  
#r=get_file(expe,date,cblock,ckind) #EXPERIENCE
r=get_file(expe,date,cblock,ckind,cfields,cmodel,cut,cterm)
#r_ref=get_file(xpref,date,cblock,ckind,cfields,cmodel,cut,cterm)

r_anaprod=get_file(expe,date,cblock,ckind,cfields,cmodel,cut_for,cterm)
r_for=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_for)
r_for_moins=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm_moins)

r_for_init=get_file(expe,date,cblock_for,ckind_for,cfields,cmodel,cut_for,cterm)

r_mer=get_file(expe,date_valid,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm)
r_mer_init=get_file(expe,date,cblock_oc,ckind_oc,cfields_oc,cmodel_oc,cut_oc,cterm)



for i in range(len(name)):
    print(i,name[i])
    
    ## LOAD FIELDS
    f=r.readfield(name[i])
    #f_ref=r_ref.readfield(name[i])
    f_mer=r_mer.readfield(name_mer[i])-273.15
    f_for=r_for.readfield(name[i])
    f_for_moins=r_for_moins.readfield(name[i])
    f_for_init=r_for_init.readfield(name[i])
    f_mer_init=r_mer_init.readfield(name_mer[i])-273.15
    f_anaprod=r_anaprod.readfield(name[i])

     
    ## LON/LAT
    X=r.geometry.get_lonlat_grid()  #lon et lat
    lon=X[0]
    lat=X[1]
    Y=r.geometry.get_lonlat_grid()
    lon_xp=Y[0]
    lat_xp=Y[1]
    
    ## MASK SUR LON/LAT
    # mask=(lon==6.24564464)&(lat==59.00403393) #masques sur lon et lat
    # mask=(lon==5.139356467619843)&(lat==59.39070739840856) #masques sur lon et lat
    #mask=(lon==5.139356467619843)&(lat==59.39070739840856)
    # mask=(lon==6.133242475917926)&(lat==59.40401122652407)
    
    # lon_zone pt 6.133242475917926
    # lat zone pt 59.40401122652407

    mask=(lon>long+0.05)&(lon<long-0.15)&(lat>lati+0.05)&(lat<lati+0.4) #masques sur lon et lat OK (5.13, 59,39) NOVEGE
    mask=(lon<long)&(lon>long-0.5)&(lat>lati)&(lat<lati+0.5) #(-2.1379, 46.9033) France
    # maskmer=(lon>long)&(lon<long+1)&(lat>lati)&(lat<lati+1)&(f_mer.data>200)&(f_mer.data<400) 
    #Long negatif
    # mask=(lon<long)&(lon>long-1)&(lat>lati)&(lat<lati+1) #masques sur lon et lat

    # #lon=np.array(lon)
    # #lat=np.array(lat)
    
    
    lon_zone=[lon[mask]]
    lat_zone=[lat[mask]]
    lon_zone=np.array(lon_zone)
    lat_zone=np.array(lat_zone)
    
    lat_pt=lat_zone[0,24]
    lon_pt=lon_zone[0,24]

    # CHAMPS DE LA ZONE
    f_total=[f.data[mask]]
    f_total=np.array(f_total)

    f_mercator=[f_mer.data[mask]]
    f_mercator=np.array(f_mercator)
    
    f_mer_total=[f_mer.data[mask]]
    f_mer_total=np.array(f_mer_total)
           
    f_for_total=[f_for.data[mask]]
    f_for_total=np.array(f_for_total)
    
    f_for_moins_total=[f_for_moins.data[mask]]
    f_for_moins_total=np.array(f_for_moins_total)
    
    f_for_init_total=[f_for_init.data[mask]]
    f_for_init_total=np.array(f_for_init_total)
    
    f_mer_init_total=[f_mer_init.data[mask]]
    f_mer_init_total=np.array(f_mer_init_total)
    
    f_prod_total=[f.data[mask]]
    f_prod_total=np.array(f_total)
    
    f_anaprod_total=[f_anaprod.data[mask]]
    f_anaprod_total=np.array(f_anaprod_total)
  

    
    # CHAMPS EN UN PT LON/LAT
    f_pt=[f_total[0,24]]
    f_pt=np.array(f_pt)
    
    ## Pareil si plusieurs pts prendre [0,2]
    ## Prend un seul point [0,0]
    f_mer_pt=[f_mer_total[0,24]]
    f_mer_pt=np.array(f_mer_pt)
  
    f_for_pt=[f_for_total[0,24]]
    f_for_pt=np.array(f_for_pt)
    
    f_for_moins_pt=[f_for_moins_total[0,24]]
    f_for_moins_pt=np.array(f_for_moins_pt)
    
    f_for_init_pt=[f_for_init_total[0,24]]
    f_for_init_pt=np.array(f_for_init_pt)

    f_mer_init_pt=[f_mer_init_total[0,24]]
    f_mer_init_pt=np.array(f_mer_init_pt)
    
    f_anaprod_pt=[f_anaprod_total[0,24]]
    f_anaprod_pt=np.array(f_anaprod_pt)

    f_prod_pt=[f_prod_total[0,24]]
    f_prod_pt=np.array(f_prod_pt)
    
     
    # # CHAMPS EN UN PT LON/LAT extract avec lon lat
    # f_pt=[f_total[0]]
    # f_pt=np.array(f_pt)
    
    # ## Pareil si plusieurs pts prendre [0,2]
    # ## Prend un seul point [0,0]
    # f_mer_pt=[f_mer_total[0]]
    # f_mer_pt=np.array(f_mer_pt)
  
    # f_for_pt=[f_for_total[0]]
    # f_for_pt=np.array(f_for_pt)
    
    # f_for_moins_pt=[f_for_moins_total[0]]
    # f_for_moins_pt=np.array(f_for_moins_pt)
    
    # f_for_init_pt=[f_for_init_total[0]]
    # f_for_init_pt=np.array(f_for_init_pt)

    # f_mer_init_pt=[f_mer_init_total[0]]
    # f_mer_init_pt=np.array(f_mer_init_pt)
    
    # f_anaprod_pt=[f_anaprod_total[0]]
    # f_anaprod_pt=np.array(f_anaprod_pt)

    # f_prod_pt=[f_prod_total[0]]
    # f_prod_pt=np.array(f_prod_pt)
    
    
    
    

    ## Append valeurs de T du profil vertical
    x=np.append(x,f_pt,axis=0)
    x_mer=np.append(x_mer,f_mer_pt,axis=0) #### Probleme ici en k
    
    x_for=np.append(x_for,f_for_pt,axis=0)
    
    x_for_moins=np.append(x_for_moins,f_for_moins_pt,axis=0)

    x_for_init=np.append(x_for_init,f_for_init_pt,axis=0)
    
    x_mer_init=np.append(x_mer_init,f_mer_init_pt,axis=0)
    
    x_anaprod=np.append(x_anaprod,f_anaprod_pt,axis=0)

    ## Fill tableau
    # T_stats=np.array([[lon_pt,lat_pt,f_pt,f_mer_pt,f_for_pt,f_for_init_pt,f_mer_init,f_anaprod]])
    # T_stats=np.array([[lon_pt,lat_pt,f_pt,f_mer_pt,f_mer_init_pt,f_for_pt,f_for_moins_pt,f_for_init]])
    # #print('T_stats',T_stats)
    # T=np.append(T,T_stats,axis=0)
    # #print('T',T)

    # # Enregistrer un tableau de donnéées format txt
    # table1 = chemin + 'Table_stats_{}.txt'.format(expe)
    # np.savetxt(table1, T, fmt='%s')


# # Profil vertical TempOC
# plt.figure()
# print('x_for uniforme',x_for)
# print('x GL55',x)
# plt.title('Prévi vs Analyse'+'-'+str(date)[0:10]+'-'+str(date)[11:13]+'_'+str(cterm_for)+'UTC'+'{}_{}'.format(long,lati))
# plt.xlabel('Temperature en °C')
# plt.ylabel('Niveaux')
# plt.plot(x, level,linewidth=2, label = str(expe)+'_ana_assim',color=color_GKOR, marker='.')
# plt.plot(x_mer, level,linewidth=2,label= 'Mercator',color=color_GN3C,linestyle='-.')
# plt.plot(x_mer_init, level,linewidth=1,label= 'Mercator init',color=color_GN3C)
# plt.plot(x_for,level,linewidth=2,label=str(expe)+'_for_P96',color=color_GKOR,linestyle='-.')
# plt.plot(x_for_init,level,linewidth=2,label=str(expe)+'_for_init_P0',color=color_GKOR,linestyle='--')
# plt.plot(x_anaprod,level,linewidth=1,label=str(expe)+'_ana_prod',color=color_GKOR,marker='+')
# #plt.plot(x_ref,level,linewidth=1,label=str(xpref)+'_ana_assim',color=color_xpref,linestyle='dotted')
# plt.legend()

# figname =  chemin+name[i]+'_previ-analyse'+str(expe)+'_'+str(date)[0:10]+'-'+str(date)[11:13]+'_'+str(cterm_for)+'UTC'+'_{}_{}.png'.format(long,lati)
# plt.savefig(figname,dpi=100, format='png')
# #print(x)
# #r.close()


# ## Profil vertical u et v
# plt.figure()
# print('x_for uniforme',x_for)
# print('x GL55',x)
# plt.title('Prévi vs Analyse'+'-'+str(date)[0:10]+'-'+str(date)[11:13]+'_'+str(cterm_for)+'UTC'+'{}_{}'.format(long,lati))
# plt.xlabel('courant u (m.s-1)')
# plt.ylabel('Niveaux')
# plt.plot(x, level, label = str(expe)+'_ana_assim',color=color_GKOR, marker='.')
# # plt.plot(x_mer, level,linewidth=2,label= 'Mercator',color=color_GN3C,linestyle='-.')
# # plt.plot(x_mer_init, level,linewidth=1,label= 'Mercator init',color=color_GN3C)
# plt.plot(x_for,level,linewidth=2,label=str(expe)+'_for_P96',color=color_GKOR,linestyle='-.')
# plt.plot(x_for_init,level,linewidth=2,label=str(expe)+'_for_init_P0',color=color_GKOR,linestyle='--')
# plt.plot(x_anaprod,level,linewidth=2,label=str(expe)+'_ana_prod',color=color_GKOR,marker='+')
# #plt.plot(x_ref,level,linewidth=1,label=str(xpref)+'_ana_assim',color=color_xpref,linestyle='dotted')

# plt.legend()



# ## Profil vertical TKE
# plt.figure()
# print('x_for uniforme',x_for)
# print('x GL55',x)
# plt.title('Prévi vs Analyse'+'-'+str(date)[0:10]+'-'+str(date)[11:13]+'_'+str(cterm_for)+'UTC'+'{}_{}'.format(long,lati))
# plt.xlabel('TKE (m2.s-2)')
# plt.ylabel('Niveaux')
# plt.plot(x, level, label = str(expe)+'_ana_assim',color=color_GKOR, marker='.')
# # plt.plot(x_mer, level,linewidth=2,label= 'Mercator',color=color_GN3C,linestyle='-.')
# # plt.plot(x_mer_init, level,linewidth=1,label= 'Mercator init',color=color_GN3C)
# plt.plot(x_for,level,linewidth=2,label=str(expe)+'_for_P96',color=color_GKOR,linestyle='-.')
# plt.plot(x_for_init,level,linewidth=2,label=str(expe)+'_for_init_P0',color=color_GKOR,linestyle='--')
# plt.plot(x_anaprod,level,linewidth=2,label=str(expe)+'_ana_prod',color=color_GKOR,marker='+')
# #plt.plot(x_ref,level,linewidth=1,label=str(xpref)+'_ana_assim',color=color_xpref,linestyle='dotted')

# plt.legend()

# Profil vertical TempOC CLEAN
plt.figure()
plt.title('Evaluation of temperature vertical profile'+' '+str(date)[0:10]+'-'+str(date)[11:13]+', lon:{},lat:{}'.format(np.round(lon_pt,3),np.round(lat_pt,3)), fontsize=13, y=1.05)
plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Depth (m)', fontsize=14)


plt.plot(x, level,linewidth=2, label = str(expe)+'_ana_assim - 09/09',color=color_expe) # ANALYSE cut off long
plt.plot(x_anaprod,level,linewidth=1,label=str(expe)+'_ana_prod - 09/09',color=color_expe,linestyle='-.')
#plt.plot(x_for_init,level,linewidth=2,label=str(expe)+'_for_init_P0',color=color_expe,linestyle='-.') # P0, Analyse cut off court

plt.plot(x_for_moins,level,linewidth=1,label=str(expe)+'_P6_prod 08/09 18UTC',color=color_expe, linestyle='dotted')
plt.plot(x_for,level,linewidth=1,label=str(expe)+'_P48_prod - 11/09',color=color_expe, linestyle='--') # P48



plt.plot(x_mer[0:12], level[0:12],linewidth=1,label= 'Mercator - 09/09',color=color_mer,linestyle='--') # Mercator a date valide
plt.plot(x_mer_init[0:12], level[0:12],linewidth=1,label= 'Mercator init - 11/09',color=color_mer) # Mercator a date init
#plt.plot(x_ref,level,linewidth=1,label=str(xpref)+'_ana_assim',color=color_xpref,linestyle='dotted')

# plt.legend(loc='upper right', bbox_to_anchor=(1.75,0.98), fontsize=12)
plt.legend(loc='best', fontsize=12)
figname =  Dir_im+name[i]+'_previ-analyse'+str(expe)+'_'+str(date)[0:10]+'-'+str(date)[11:13]+'_'+str(cterm_for)+'UTC'+'_{}_{}.png'.format(np.round(lon_pt,3),np.round(lat_pt,3))
plt.savefig(figname,dpi=300, format='png')

# # Profil vertical TempOC
# plt.figure()
# print('x_for uniforme',x_for)
# print('x GL55',x)
# plt.title('Prévi vs Analyse'+'-'+str(date)[0:10]+'-'+str(date)[11:13]+'_'+str(cterm_for)+'UTC'+'{}_{}'.format(long,lati))
# plt.xlabel('Salinité')
# plt.ylabel('Niveaux')
# plt.plot(x, level,linewidth=2, label = str(expe)+'_ana_assim',color=color_GM6E, marker='.')
# plt.plot(x_mer, level,linewidth=2,label= 'Mercator',color=color_GN3C,linestyle='-.')
# plt.plot(x_mer_init, level,linewidth=1,label= 'Mercator init',color=color_GN3C)
# plt.plot(x_for,level,linewidth=2,label=str(expe)+'_for_P96',color=color_GM6E,linestyle='-.')
# plt.plot(x_for_init,level,linewidth=2,label=str(expe)+'_for_init_P0',color=color_GM6E,linestyle='--')
# plt.plot(x_anaprod,level,linewidth=1,label=str(expe)+'_ana_prod',color=color_GM6E,marker='+')
# #plt.plot(x_ref,level,linewidth=1,label=str(xpref)+'_ana_assim',color=color_xpref,linestyle='dotted')
# plt.legend()