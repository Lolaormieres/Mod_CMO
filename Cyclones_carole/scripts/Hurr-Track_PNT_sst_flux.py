#!/usr/bin/env python3
# -*- coding: utf-8 -*-



### Import

import epygram
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os 
import csv
from datetime import datetime,timedelta
from vortex import toolbox
import common, olive
import common.util.usepygram
import usevortex


epygram.init_env()


from fonctions2 import nombre_TC, recup_data_ibtrack, latlon_param_TC_parmembre_pour1ech#, recup_size_ibtrack_parTC

#####################################################################
# Initializations
#####################################################################
#pathrun='tmpdir'

def get_file(xpid,date,term):
    #MODEL
     geometry=domaine
     fp = {'experiment':xpid,
          'kind':'gridpoint',
          'block':'forecast',
          'date':date,
          'term':term,
          'geometry':geometry,#vortex.data.geometries.get(tag='global1798'),
          'format':'grib',
#          'local':'tmp2.fic',
          'shouldfly': True,
          'cutoff':'production',
          'vapp':'arpege',
          'vconf':'4dvarfr',
          'model':'arpege',
          'namespace':'vortex.multi.fr',
          'origin':'historic',
          'nativefmt':'grib',
          'now': True,
          }
     fcref = toolbox.input(**fp)
     r = fcref[0].contents.data
     #return {int(r.provider.member): r.contents.data for r in fcref}
     return r

# Récupération des variables dans la namelist
fd=open('Hurr-Track.nam','r')
data_nam=np.loadtxt('Hurr-Track.nam',dtype='str')[:]


annee=data_nam[0]
mois=data_nam[1]
jour=data_nam[2]
heure=data_nam[3]
pdt=int(data_nam[4])
model=data_nam[5]
member=data_nam[6]
domaine=data_nam[7]
ech_min=int(data_nam[8])
ech_max=int(data_nam[9])
expe=data_nam[10]

print('****************',annee, mois, jour, heure,pdt)

date=datetime(int(annee),int(mois),int(jour),int(heure))
print(date)

print(member)
  
if(domaine=='glob025'):
  zdeltax = 25000.
  zdeltay = 25000.
  iie=1440
  ije=721


if(domaine=='glob01'):
  zdeltax = 10000.
  zdeltay = 10000.
  iie=3600
  ije=1801

"""  
if(domaine=='glob05'):
  zdeltax = 50000.
  zdeltay = 50000.
  iie=720
  ije=361  
"""

list_ech=np.arange(ech_min,ech_max+pdt,pdt)
ech=len(list_ech)
list_tc=[]

speedguess = 50.
xradguess = 600. # en km
xradguessVmax = 300. # en km
nphil=33  # number of azimuts
xboxwind=30. # 
zspeedmin=100.
iib=0
ijb=0

zdeltar = max(zdeltax,zdeltay)
zdphi = np.math.pi / 16.

iradmax0 = int((xradguess*1000. / zdeltar) + 2)
iradVmax0 = int((xradguessVmax*1000. / zdeltar) + 2)

liste_bassin=["NA", "EP", "WP", "SP", "SI", "NI"]

g=open("nbcyclone.txt",'r')
text=csv.reader(g, delimiter=',')
tab_cycl=list(text)
g.close()
########################################
# Lancement
########################################


for jc in range(ech):

  print('On traite l\'échéance '+str(list_ech[jc]).zfill(2))
  #lecture 
  path = os.getcwd()
  print("Le répertoire courant est : " + path)
  f=get_file(expe,date,int(list_ech[jc]))
#        f = epygram.formats.resource('grid.'+model+'-forecast.'+domaine+'+'+str(list_ech[jc]).zfill(4)+':00.grib', 'r', 'GRIB')
  path = os.getcwd()
  print("Le répertoire courant est : " + path)

# lecture des variables
  mslp = f.readfield({'name': 'Pressure reduced to MSL'})
  mslp.operation('/',100.)

  u10=f.readfield({'name':'10 metre U wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
  v10=f.readfield({'name':'10 metre V wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
  u850=f.readfield({'name':'U component of wind', 'typeOfFirstFixedSurface':100, 'level': 850})
  v850=f.readfield({'name':'V component of wind', 'typeOfFirstFixedSurface':100, 'level': 850})
  u925=f.readfield({'name':'U component of wind', 'typeOfFirstFixedSurface':100, 'level': 925})
  v925=f.readfield({'name':'V component of wind', 'typeOfFirstFixedSurface':100, 'level': 925})
  u50=f.readfield({'name':'U component of wind', 'typeOfFirstFixedSurface':103, 'level': 50})
  v50=f.readfield({'name':'V component of wind', 'typeOfFirstFixedSurface':103, 'level': 50})
  u100=f.readfield({'name':'100 metre U wind component', 'typeOfFirstFixedSurface':103, 'level': 100})
  v100=f.readfield({'name':'100 metre V wind component', 'typeOfFirstFixedSurface':103, 'level': 100})
  u500=f.readfield({'name':'U component of wind', 'typeOfFirstFixedSurface':103, 'level': 500})
  v500=f.readfield({'name':'V component of wind', 'typeOfFirstFixedSurface':103, 'level': 500})
  ugust1h=f.readfield({'name':'10 metre eastward wind gust since previous post-processing', 'stepRange':str(list_ech[jc]-1)+'-'+str(list_ech[jc])})
  vgust1h=f.readfield({'name':'10 metre northward wind gust since previous post-processing', 'stepRange':str(list_ech[jc]-1)+'-'+str(list_ech[jc])})
  ugust3h=f.readfield({'name':'10 metre eastward wind gust since previous post-processing', 'stepRange':str(list_ech[jc]-3)+'-'+str(list_ech[jc])})
  vgust3h=f.readfield({'name':'10 metre northward wind gust since previous post-processing', 'stepRange':str(list_ech[jc]-3)+'-'+str(list_ech[jc])})
  sstr=f.readfield({'shortName':'sot'})
  flatr=f.readfield({'shortName':'slhf'})
  fsenr=f.readfield({'shortName':'sshf'})
  if (jc>0):
      fb=get_file(expe,date,int(list_ech[jc-1]))
      flatrb=fb.readfield({'shortName':'slhf'})
      fsenrb=fb.readfield({'shortName':'sshf'})

  for ligne in range(len(tab_cycl)) :

      nbr_tc=int(tab_cycl[ligne][0])
      ibbassin=tab_cycl[ligne][1]
    
#### Loop on the number of analysed tropical cyclone

      for jname in range(nbr_tc):

        zradius=np.zeros(iradmax0, dtype=int)
        ztime=np.zeros(ech, dtype=object)
        ztime_basis=np.zeros(1,dtype=object)
        zspeedmin=100.

      # vents 850,925 hPa
        zmod850cyl = np.zeros((ech,nphil,iradmax0), dtype=float)
        zmod850cylMAX = np.zeros((ech,iradmax0), dtype=float)
    
        zmod925cyl = np.zeros((ech,nphil,iradmax0), dtype=float)
        zmod925cylMAX = np.zeros((ech,iradmax0), dtype=float)

      # vents 10,50,100,500 m
        zmod10mcyl = np.zeros((ech,nphil,iradmax0), dtype=float)
        zmod10mcylNE = np.zeros((ech,iradmax0), dtype=float)
        zmod10mcylNW = np.zeros((ech,iradmax0), dtype=float)
        zmod10mcylSW = np.zeros((ech,iradmax0), dtype=float)
        zmod10mcylSE = np.zeros((ech,iradmax0), dtype=float)
        zmod10mcylMAX = np.zeros((ech,iradmax0), dtype=float)
 
        zmod50mcyl = np.zeros((ech,nphil,iradmax0), dtype=float)
        zmod50mcylMAX = np.zeros((ech,iradmax0), dtype=float)

        zmod100mcyl = np.zeros((ech,nphil,iradmax0), dtype=float)
        zmod100mcylMAX = np.zeros((ech,iradmax0), dtype=float)

        zmod500mcyl = np.zeros((ech,nphil,iradmax0), dtype=float)
        zmod500mcylMAX = np.zeros((ech,iradmax0), dtype=float)

     # rafales sur 1h, 3h avant la prévision
        zmodgust1hcyl = np.zeros((ech,nphil,iradmax0), dtype=float)
        zmodgust1hcylMAX = np.zeros((ech,iradmax0), dtype=float)

        zmodgust3hcyl = np.zeros((ech,nphil,iradmax0), dtype=float)
        zmodgust3hcylMAX = np.zeros((ech,iradmax0), dtype=float)

      # récup datas
#      tcname=list_tc[jname]
        j2name=2+jname
        tcname=tab_cycl[ligne][j2name]
        data_guess=recup_data_ibtrack(annee,mois, jour, heure, ibbassin, tcname)
    
        latguess=np.asarray(data_guess[:,8],dtype='float' )
        longuess=np.asarray(data_guess[:,9],dtype='float')
    
        if(ibbassin=='NA' or ibbassin=='EP'):
          print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % latguess)+' °N et '+str("% 4.2f" % (longuess*-1.))+' °W')
        if(ibbassin=='SI' or ibbassin=='SP'):
          print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % (latguess*-1.))+' °S et '+str("% 4.2f" % longuess)+' °E')       
        if(ibbassin=='WP' or ibbassin=='NI'):
          print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % latguess)+' °N et '+str("% 4.2f" % (longuess))+' °E')

# boucle sur les échéances
        ztime_basis[0]=mslp.validity.getbasis()
        ztime_term=mslp.validity.term()
        ztime[jc]=mslp.validity[0].get()

        
        ##########################################################################
        # Détermination de la zone de caccul des conditions cycloniques 
        # (on considère un cyclone potentiel se déplaçant à une vitesse speedguess
        ##########################################################################

        if(jc-1>=0):

           data0 = np.loadtxt(tcname+'_'+annee+mois+jour+'_'+heure+'.csv',dtype="str", delimiter=',')
           print("taille de data0",np.size(data0))
           print(data0)
           if np.size(data0) == 25 :
             latpos_1 = np.asarray(data0[3], dtype='float')
             lonpos_1 = np.asarray(data0[4], dtype='float')
           else :
             data = data0[(data0[:,1] == member) & (data0[:,2] == str(list_ech[jc-1])), :]
             latpos_1 = np.asarray(data[:,3], dtype='float')
             lonpos_1 = np.asarray(data[:,4], dtype='float')

        else:
           latpos_1 = latguess
           lonpos_1 = longuess

        distguess = (pdt * speedguess) / 111.
        #
        latmaxpos = latpos_1 + distguess
        latminpos = latpos_1 - distguess
        lonmaxpos = lonpos_1 + distguess
        lonminpos = lonpos_1 - distguess
        #
        niboxinf = int(mslp.geometry.ll2ij(lonminpos,latminpos)[0])
        niboxsup = int(mslp.geometry.ll2ij(lonmaxpos,latminpos)[0])
        njboxinf = int(mslp.geometry.ll2ij(lonminpos,latminpos)[1])
        njboxsup = int(mslp.geometry.ll2ij(lonminpos,latmaxpos)[1])


        #####################################################################
        # Détermination du centre cyclonique iicen et ijcen 
        # (minimum de pression réduite au niveau de la mer)
        #####################################################################

        mslpmin=1050
        for j in range(njboxinf,njboxsup):
            for i in range(niboxinf,niboxsup):
                if mslp.data[j,i]<mslpmin:
                   mslpmin=mslp.data[j,i]
                   iicen = i
                   ijcen = j
#

        print('mslpmin1=',mslpmin)
#        print('njboxinf:njboxsup,niboxinf:niboxsup=',njboxinf,njboxsup,niboxinf,niboxsup)
#        print('y a quoi la dedans=',mslp.data[njboxinf:njboxsup,niboxinf:niboxsup])
        mslpmin=np.min(mslp.data[njboxinf:njboxsup,niboxinf:niboxsup])
        print('mslpmin2=',mslpmin)
#        print('iicen,ijcen attendus=',np.unravel_index(np.argmin(mslp.data[njboxinf:njboxsup,niboxinf:niboxsup], axis=None), mslp.data[njboxinf:njboxsup,niboxinf:niboxsup].shape))
#        print('iicen,ijcen trouvés=',iicen,ijcen)
#        iicen,ijcen=np.unravel_index(np.argmin(mslp.data[njboxinf:njboxsup,niboxinf:niboxsup], axis=None), mslp.data[njboxinf:njboxsup,niboxinf:niboxsup].shape)
#        iicen,ijcen=np.argmin(mslp.data[njboxinf:njboxsup,niboxinf:niboxsup],axis=0)
        centre_ll=mslp.geometry.ij2ll(iicen,ijcen)
        print('centre_ll=',centre_ll)
#   
        #####################################################################
        # Détermination du centre cyclonique iicen et ijcen 
        # (minimum de vent dans l'oeil)
        #####################################################################
#


        mod10m = np.sqrt(u10.data[:,:]**2 + v10.data[:,:]**2)
#
        iavgwindi = int(round(xboxwind*1000./zdeltax))
        iavgwindj = int(round(xboxwind*1000./zdeltay))
#
        for jj in range(max(ijcen-iavgwindj,1),min(ijcen+iavgwindj,ije)):
            for ji in range(max(iicen-iavgwindi,1),min(iicen+iavgwindi,iie)):
                zspeed=np.sqrt(u10.data[jj,ji]**2. + v10.data[jj,ji]**2.)
                if zspeed<= zspeedmin:
                   zspeedmin=zspeed
                   iicen2 = ji
                   ijcen2 = jj
#  
        iicen=iicen2
        ijcen=ijcen2
#
        centre_ll=mslp.geometry.ij2ll(iicen,ijcen)

        if(ibbassin=='NA' or ibbassin=='EP'):
           print('A l\'échéance '+str(list_ech[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % centre_ll[1])+' °N et'+str("% 5.2f" % ((centre_ll[0]-360.)*-1)+' °W'))
        if(ibbassin=='SI' or ibbassin=='SP'):
           print('A l\'échéance '+str(list_ech[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % (centre_ll[1]*-1.))+' °S et'+str("% 5.2f" % centre_ll[0]+' °E'))   
        if(ibbassin=='WP' or ibbassin=='NI'):
           print('A l\'échéance '+str(list_ech[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % centre_ll[1])+' °N et'+str("% 5.2f" % centre_ll[0]+' °E'))  


        # sst et flux lat
        sst=sstr.data[ijcen,iicen]
        if (jc>0):
            flat=(flatr.data[ijcen,iicen]-flatrb.data[ijcen,iicen])/(pdt*3600)
            fsen=(fsenr.data[ijcen,iicen]-fsenrb.data[ijcen,iicen])/(pdt*3600)
        else:
            flat=flatr.data[ijcen,iicen]/(pdt*3600)
            fsen=fsenr.data[ijcen,iicen]/(pdt*3600)

        ##############################################################################
        # Ajout des différents paramètres de vents pour le calcul du vent max 
        # vents 850,925 hPa ; 50,100,500m ; rafales sur 1, 3h précédent la prévision
        ##############################################################################

#
        
        mod850 = np.sqrt(u850.data[:,:]**2 + v850.data[:,:]**2)

#
        
        mod925 = np.sqrt(u925.data[:,:]**2 + v925.data[:,:]**2)

#

        mod50m = np.sqrt(u50.data[:,:]**2 + v50.data[:,:]**2)
        
#

        mod100m = np.sqrt(u100.data[:,:]**2 + v100.data[:,:]**2)

#

        mod500m = np.sqrt(u500.data[:,:]**2 + v500.data[:,:]**2)

#
#        ugust1h=f.readfield({'name':'u-component of wind (gust)', 'stepRange':str(list_ech[jc]-1)+'-'+str(list_ech[jc])})
#        vgust1h=f.readfield({'name':'v-component of wind (gust)', 'stepRange':str(list_ech[jc]-1)+'-'+str(list_ech[jc])})

        modgust1h = np.sqrt(ugust1h.data[:,:]**2 + vgust1h.data[:,:]**2) 
     
#
#        ugust3h=f.readfield({'name':'u-component of wind (gust)', 'stepRange':str(list_ech[jc]-3)+'-'+str(list_ech[jc])})
#        vgust3h=f.readfield({'name':'v-component of wind (gust)', 'stepRange':str(list_ech[jc]-3)+'-'+str(list_ech[jc])})

        modgust3h = np.sqrt(ugust3h.data[:,:]**2 + vgust3h.data[:,:]**2)  


#
        #####################################################################
        # Calcul des paramètres en géométrie cylindrique
        # (R=1 corresponds au centre du cyclone)
        #####################################################################

#           
        for jr in range(1,iradmax0):
            zradius[jr] = jr * zdeltar
# 
        zxi0=mslp.geometry.ij2xy(iicen,ijcen)[0]*100000. 
        zyj0=mslp.geometry.ij2xy(iicen,ijcen)[1]*100000. 
        zx00=mslp.geometry.ij2xy(1,ijcen)[0]*100000. 
        zy00=mslp.geometry.ij2xy(iicen,1)[1]*100000. 
# 
        #initialisations pour les calculs d extensions radiales de vent
        zrmwmaxNE=0.
        zrmwmaxNW=0.
        zrmwmaxSW=0.
        zrmwmaxSE=0.
        irmwNE=0
        irmwNW=0
        irmwSW=0
        irmwSE=0

        for jr in range(iradmax0):
            for jphi in range(nphil):
                zphi = (jphi) * zdphi
                zxk  = zradius[jr] * np.cos(zphi) + zxi0
                zyk  = zradius[jr] * np.sin(zphi) + zyj0
                iix  = int(((zxk-zx00) / zdeltax) + 1 )
                iiy  = int(((zyk-zy00) / zdeltay) + 1)
#
                if (iix>iib) and ((iix+1)<iie) and (iiy>ijb) and ((iiy+1)<ije): 
                   zxk = (zxk-mslp.geometry.ij2xy(iix,iiy)[0]*100000.) / zdeltax 
                   zyk = (zyk-mslp.geometry.ij2xy(iix,iiy)[1]*100000.) / zdeltay 
                                  
                   zmod850cyl[jc,jphi,jr]=mod850.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                     mod850.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                     mod850.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                     mod850.data[iiy+1,iix+1]*zxk*zyk

                   zmod925cyl[jc,jphi,jr]=mod925.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                     mod925.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                     mod925.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                     mod925.data[iiy+1,iix+1]*zxk*zyk

                                     
                   zmod10mcyl[jc,jphi,jr]=mod10m.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                     mod10m.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                     mod10m.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                     mod10m.data[iiy+1,iix+1]*zxk*zyk  

                   zmod50mcyl[jc,jphi,jr]=mod50m.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                     mod50m.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                     mod50m.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                     mod50m.data[iiy+1,iix+1]*zxk*zyk                  

                   zmod100mcyl[jc,jphi,jr]=mod100m.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                     mod100m.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                     mod100m.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                     mod100m.data[iiy+1,iix+1]*zxk*zyk                  

                   zmod500mcyl[jc,jphi,jr]=mod500m.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                     mod500m.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                     mod500m.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                     mod500m.data[iiy+1,iix+1]*zxk*zyk     


                   zmodgust1hcyl[jc,jphi,jr]=modgust1h.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                     modgust1h.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                     modgust1h.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                     modgust1h.data[iiy+1,iix+1]*zxk*zyk          

                   zmodgust3hcyl[jc,jphi,jr]=modgust3h.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                     modgust3h.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                     modgust3h.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                     modgust3h.data[iiy+1,iix+1]*zxk*zyk         
             
                

                   

                   zmod850cyl[jc,jphi,jr] = zmod850cyl[jc,jphi,jr] * 1.944   # Conversion en kts
                   zmod925cyl[jc,jphi,jr] = zmod925cyl[jc,jphi,jr] * 1.944   # Conversion en kts
                   zmod10mcyl[jc,jphi,jr] = zmod10mcyl[jc,jphi,jr] * 1.944   # Conversion en kts
                   zmod50mcyl[jc,jphi,jr] = zmod50mcyl[jc,jphi,jr] * 1.944   # Conversion en kts
                   zmod100mcyl[jc,jphi,jr] = zmod100mcyl[jc,jphi,jr] * 1.944   # Conversion en kts
                   zmod500mcyl[jc,jphi,jr] = zmod500mcyl[jc,jphi,jr] * 1.944   # Conversion en kts
                   zmodgust1hcyl[jc,jphi,jr] = zmodgust1hcyl[jc,jphi,jr] * 1.944   # Conversion en kts
                   zmodgust3hcyl[jc,jphi,jr] = zmodgust3hcyl[jc,jphi,jr] * 1.944   # Conversion en kts
                   
                   #
                else:
            
                   zmod850cyl[jc,jphi, jr]=-999.
                   zmod925cyl[jc,jphi, jr]=-999.
         
                   zmod10mcyl[jc,jphi, jr]=-999.
                   zmod50mcyl[jc,jphi, jr]=-999.
                   zmod100mcyl[jc,jphi, jr]=-999.
                   zmod500mcyl[jc,jphi, jr]=-999.

                   zmodgust1hcyl[jc,jphi, jr]=-999.
                   zmodgust3hcyl[jc,jphi, jr]=-999.


            zmod850cylMAX[jc,jr] = np.max(zmod850cyl[jc,:,jr])
            zmod925cylMAX[jc,jr] = np.max(zmod925cyl[jc,:,jr])
            
            zmod10mcylMAX[jc,jr] = np.max(zmod10mcyl[jc,:,jr])
            zmod50mcylMAX[jc,jr] = np.max(zmod50mcyl[jc,:,jr])
            zmod100mcylMAX[jc,jr] = np.max(zmod100mcyl[jc,:,jr])
            zmod500mcylMAX[jc,jr] = np.max(zmod500mcyl[jc,:,jr])

            zmodgust1hcylMAX[jc,jr] = np.max(zmodgust1hcyl[jc,:,jr])
            zmodgust3hcylMAX[jc,jr] = np.max(zmodgust3hcyl[jc,:,jr])
   
# 
        ###################################################################
        #Calcul les extensions radiales de vent des cyclones dans la PEARP.
        ###################################################################
        
            icount=0
            for jphi in range(0,9):
                if(zmod10mcyl[jc,jphi,jr]!=-999.):
                    zmod10mcylNE[jc,jr] = zmod10mcylNE[jc,jr] + zmod10mcyl[jc,jphi,jr]
                    icount = icount + 1
                else:
                    zmod10mcylNE[jc,jr]=-999.
            if icount == len(range(0,9)):
                zmod10mcylNE[jc,jr] = zmod10mcylNE[jc,jr] / icount
            else:
                zmod10mcylNE[jc,jr] = -999.

                # range(8,17) pour zmod10mcylNW
            icount=0
            for jphi in range(8,17):
                if(zmod10mcyl[jc,jphi,jr]!=-999.):
                    zmod10mcylNW[jc,jr] = zmod10mcylNW[jc,jr] + zmod10mcyl[jc,jphi,jr]
                    icount = icount + 1
                else:
                    zmod10mcylNW[jc,jr]=-999.
            if icount == len(range(8,17)):
                zmod10mcylNW[jc,jr] = zmod10mcylNW[jc,jr] / icount
            else:
                zmod10mcylNW[jc,jr] = -999.

                # range(16,25) pour zmod10mcySW
            icount=0
            for jphi in range(16,25):
                if(zmod10mcyl[jc,jphi,jr]!=-999.):
                    zmod10mcylSW[jc,jr] = zmod10mcylSW[jc,jr] + zmod10mcyl[jc,jphi,jr]
                    icount = icount + 1
                else:
                    zmod10mcylSW[jc,jr]=-999.
            if icount == len(range(16,25)):
                zmod10mcylSW[jc,jr] = zmod10mcylSW[jc,jr] / icount
            else:
                zmod10mcylSW[jc,jr] = -999.

                # range(24,33) pour zutcylSE
            icount=0
            for jphi in range(24,33):
                if(zmod10mcyl[jc,jphi,jr]!=-999.):
                    zmod10mcylSE[jc,jr] = zmod10mcylSE[jc,jr] + zmod10mcyl[jc,jphi,jr]
                    icount = icount + 1
                else:
                    zmod10mcylSE[jc,jr]=-999.
            if icount == len(range(24,33)):
                zmod10mcylSE[jc,jr] = zmod10mcylSE[jc,jr] / icount
            else:
                zmod10mcylSE[jc,jr] = -999.

                #
            if zmod10mcylNE[jc,jr]>zrmwmaxNE:
                zrmwmaxNE=zmod10mcylNE[jc,jr]
                irmwNE=jr

                #
            if zmod10mcylNW[jc,jr]>zrmwmaxNW:
                zrmwmaxNW=zmod10mcylNW[jc,jr]
                irmwNW=jr
                #
            if zmod10mcylSW[jc,jr]>zrmwmaxSW:
                zrmwmaxSW=zmod10mcylSW[jc,jr]
                irmwSW=jr
                #
            if zmod10mcylSE[jc,jr]>zrmwmaxSE:
                zrmwmaxSE=zmod10mcylSE[jc,jr]
                irmwSE=jr

                #Pour le seuil 34 
        zrmw34NE=-999.
        zrmw34NW=-999.
        zrmw34SW=-999.
        zrmw34SE=-999.
                #
        if (irmwNE>0):
            if (zmod10mcylNE[jc,irmwNE]>34):
                for jr in range(irmwNE,iradmax0):
                    if (zmod10mcylNE[jc,jr]<34 and zmod10mcylNE[jc,jr]!=-999.):
                        irmw34=jr
                        zrmw34NE=irmwNE*zdeltar / 1000. + (irmw34-irmwNE)*zdeltar / 1000.
                        break

        if (irmwNW>0):
            if (zmod10mcylNW[jc,irmwNW]>34):
                for jr in range(irmwNW,iradmax0):
                    if (zmod10mcylNW[jc,jr]<34 and zmod10mcylNW[jc,jr]!=-999.):
                        irmw34=jr
                        zrmw34NW=irmwNW*zdeltar / 1000. + (irmw34-irmwNW)*zdeltar / 1000.
                        break

        if (irmwSW>0):
            if (zmod10mcylSW[jc,irmwSW]>34):
                for jr in range(irmwSW,iradmax0):
                    if (zmod10mcylSW[jc,jr]<34 and zmod10mcylSW[jc,jr]!=-999.):
                        irmw34=jr
                        zrmw34SW=irmwSW*zdeltar / 1000. + (irmw34-irmwSW)*zdeltar / 1000.
                        break

        if (irmwSE>0):
            if (zmod10mcylSE[jc,irmwSE]>34):
                for jr in range(irmwSE,iradmax0):
                    if (zmod10mcylSE[jc,jr]<34 and zmod10mcylSE[jc,jr]!=-999.):
                        irmw34=jr
                        zrmw34SE=irmwSE*zdeltar / 1000. + (irmw34-irmwSE)*zdeltar / 1000.
                        break
             #Calcul de la moyenne des distances pour le seuil 34
        #if(zrmw34NE>0 and zrmw34NW>0 and zrmw34SW>0 and zrmw34SE>0):  
            #zrmw34MEAN=(zrmw34NE+zrmw34NW+zrmw34SW+zrmw34SE) / 4.
        #else : 
            #zrmw34MEAN=-999."

                #
                #Pour le seuil 64
        zrmw64NE=-999.
        zrmw64NW=-999.
        zrmw64SW=-999.
        zrmw64SE=-999.


                #
        if (irmwNE>0):
            if (zmod10mcylNE[jc,irmwNE]>64):
                for jr in range(irmwNE,iradmax0):
                    if (zmod10mcylNE[jc,jr]<64 and zmod10mcylNE[jc,jr]!=-999.):
                        irmw64=jr
                        zrmw64NE=irmw64*zdeltar / 1000.
                        break

        if (irmwNW>0):
            if (zmod10mcylNW[jc,irmwNW]>64):
                for jr in range(irmwNW,iradmax0):
                    if (zmod10mcylNW[jc,jr]<64 and zmod10mcylNW[jc,jr]!=-999.):
                        irmw64=jr
                        zrmw64NW=irmw64*zdeltar / 1000.
                        break

        if (irmwSW>0):
            if (zmod10mcylSW[jc,irmwSW]>64):
                for jr in range(irmwSW,iradmax0):
                    if (zmod10mcylSW[jc,jr]<64 and zmod10mcylSW[jc,jr]!=-999.):
                        irmw64=jr
                        zrmw64SW=irmw64*zdeltar / 1000.
                        break

        if (irmwSE>0):
            if (zmod10mcylSE[jc,irmwSE]>64):
                for jr in range(irmwSE,iradmax0):
                    if (zmod10mcylSE[jc,jr]<64 and zmod10mcylSE[jc,jr]!=-999.):
                        irmw64=jr
                        zrmw64SE=irmw64*zdeltar / 1000.
                        break
               #
        #Calcul de la moyenne des distances pour le seuil 64
        #if(zrmw64NE>0 and zrmw64NW>0 and zrmw64SW>0 and zrmw64SE>0):
            #zrmw64MEAN=(zrmw64NE+zrmw64NW+zrmw64SW+zrmw64SE) / 4.
        #else :
            #zrmw64MEAN=-999.

#        
        ###############################################################
        # Création fichiers de paramètres
        ###############################################################

        f=open(tcname+'_'+annee+mois+jour+'_'+heure+'.csv','a')

	
        if(ibbassin=='NA' or ibbassin=='EP'):
           lat=centre_ll[1]
           lon=centre_ll[0]-360.
        if(ibbassin=='SI' or ibbassin=='WP' or ibbassin=='SP' or ibbassin=='NI'):   
           lat=centre_ll[1]
           lon=centre_ll[0]
       


        # vents max
        vmax850=np.max(zmod850cylMAX[jc,:])
        vmax925=np.max(zmod850cylMAX[jc,:])

        vmax10m=np.max(zmod10mcylMAX[jc,:])
        vmax50m=np.max(zmod50mcylMAX[jc,:])
        vmax100m=np.max(zmod100mcylMAX[jc,:])
        vmax500m=np.max(zmod500mcylMAX[jc,:])

        vmaxgust1h=np.max(zmodgust1hcylMAX[jc,:])
        vmaxgust3h=np.max(zmodgust3hcylMAX[jc,:])
      
        # pression minimale
        pmin=mslpmin
        print('sst*******************************************************************',sst)


        f.write(ibbassin+','+member+','+str(list_ech[jc])+','+str("%4.2f" % lat)+','+str("%4.2f" % lon)+','+str("%2i" % vmax850)+','+str("%2i" % vmax925)+','+str("%2i" % vmax10m)+','+str("%2i" % vmax50m)+','+str("%2i" % vmax100m)+','+str("%2i" % vmax500m)+','+str("%2i" % vmaxgust1h)+','+str("%2i" % vmaxgust3h)+','+str("%4i" % pmin)+','+str("%4i" % zrmw34NE)+','+str("%4i" % zrmw34NW)+','+str("%4i" % zrmw34SW)+','+str("%4i" % zrmw34SE)+','+str("%4i" % zrmw64NE)+','+str("%4i" % zrmw64NW)+','+str("%4i" % zrmw64SW)+','+str("%4i" % zrmw64SE)+','+str("%4i" % sst)+','+str("%4i" % flat)+','+str("%4i" % fsen)+"\n")
        #f.write(ibbassin+','+member+','+str(list_ech[jc])+','+str("%4.2f" % lat)+','+str("%4.2f" % lon)+','+str("%2i" % vmax850)+','+str("%2i" % vmax925)+','+str("%2i" % vmax10m)+','+str("%2i" % vmax50m)+','+str("%2i" % vmax100m)+','+str("%2i" % vmax500m)+','+str("%2i" % vmaxgust1h)+','+str("%2i" % vmaxgust3h)+','+str("%4i" % pmin)+','+str("%4i" % zrmw34NE)+','+str("%4i" % zrmw34NW)+','+str("%4i" % zrmw34SW)+','+str("%4i" % zrmw34SE)+','+str("%4i" % zrmw34MEAN)+','+str("%4i" % zrmw64NE)+','+str("%4i" % zrmw64NW)+','+str("%4i" % zrmw64SW)+','+str("%4i" % zrmw64SE)+','+str("%4i" % zrmw64MEAN)+"\n")

        f.close()
            
  fileaeffacer='grid.arpege-forecast.glob01+'+str(list_ech[jc]).zfill(4)+':00.grib'
  effacecadir='/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr/'+expe+'/'+annee+mois+jour+'T'+heure+'00P/forecast/'
  effaceca = os.path.join(effacecadir, fileaeffacer)
  print(os.listdir(path=effacecadir))
  os.remove(effaceca)
#  if (jc>0):
#    fileaeffacer2='grid.arpege-forecast.glob01+'+str(list_ech[jc-1]).zfill(4)+':00.grib'
#    effaceca2 = os.path.join(effacecadir, fileaeffacer2)
#    os.remove(effaceca2)




