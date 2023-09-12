##!/usr/bin/env python
#!/bin/env python

from Magics.macro import *
import epygram
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from PIL import Image
from compute_TC_features_tools import polar_calc, polar_calcwind, mean_azim, read_BT, fix_TC_loc, write_model_data, tracker_run, compute_RMW, compute_Vtmax_azim, compute_windshear
from plot_TC_features_tools import locate_tc_wind_2D, locate_tc_rain_2D, locate_tc_isp_2D

epygram.init_env()




#####################################################################
# Initializations
#####################################################################


fd=open('Hurr-Track.nam','r')
data_nam=numpy.loadtxt('Hurr-Track.nam',dtype='str')[:]

annee=data_nam[0]
mois=data_nam[1]
jour=data_nam[2]
heure=data_nam[3]
model=data_nam[4]
domaine=data_nam[5]
area=data_nam[6]
bassin=data_nam[7]
tcname=data_nam[8]
echmax=int(data_nam[9])

LGRIB2=True   # Ne pas changer car c'est du grib2 depuis 2019

zdeltax = 2500. # Résolution horizontale grille Arome en m
zdeltay = 2500.


print("On se trouve sur le bassin=", bassin)


if(domaine=='CARAIB0025'):
  Hemisphere = 'N'  
  iie=945
  ije=529

if(domaine=='atllarge0025'):
  Hemisphere = 'N'  
  iie=1201
  ije=841  

if(domaine=='INDIEN0025'):
  Hemisphere = 'S'
  iie=1395 # Nbre de points de la grille de sortie 
  ije=747

if(domaine=='NCALED0025'):
  iie=521
  ije=491  



speedguess = 50. # Vitesse d'advection de la boîte pour la recherche du centre cyclonique
xradguess = 500. # Distance depuis le centre du cyclone en coordonnées polaires en km
zdeltar = max(zdeltax,zdeltay)    # Résolution radiale en polaire (même que la grille de sortie)
zdphi = numpy.math.pi / 16. # Résolution azimutale
nphil = int((2. * numpy.math.pi / zdphi)) + 1 # Nbre d'azimuts en fonction de la résolution
xboxwind=30. # Boite autour de la Pmin pour trouver le min de vent comme centre (attention désactivé pour le moment dans tracker_run)
iib=0
ijb=0

KMX=400.   # Distance souhaitée autour du centre cyclonique dans le repère cartésien (en km)
KMY=400.

NX=int(KMX/(zdeltax/1000.))    
NY=int(KMY/(zdeltay/1000.))

#   Fixer les échéances et pas de temps ici
echdeb=0
echfin=echmax+3
echstep=3
stepboxmove = echstep



ficobs='BT_HISTORY_'+tcname
#ficEPS='PE'+area+'-'+annee+mois+jour+heure+'-'+tcname

list_mb=[]
for m in range(0,17,1):
    list_mb.append(m) # Les niveaux lus dans les gribs Arome
mb=len(list_mb)


list_ech=[]
for t in range(echdeb,echfin,echstep):
    list_ech.append(t) # Les niveaux lus dans les gribs Arome
ech=len(list_ech)

list_tc=[]

##### Needed for computing wind shear
xradguessS = 500. # en km
zdeltarS = 100000. # On decoupe par anneau de 100 km
zdphiS = numpy.math.pi / 36. # On decoupe par azimut de 5 degres
iradmax0S = numpy.int64((xradguessS*1000. / zdeltarS) + 1)
nphilS=int((2. * numpy.math.pi / zdphiS)) + 1  # Nbre d'azimuts en fonction de la résolution

iradmax0 = numpy.int64((xradguess*1000. / zdeltar) + 2)



### For Visualisation
zoom= NX * zdeltax / 1000.


#####################################################################
# Initializations
#####################################################################


nbr_tc = 1   # On force la loop a 1 TC

print('Il y a '+str(nbr_tc)+' cyclones analysés pour le '+jour+'/'+mois+'/'+annee+' à '+heure+' UTC')


### Read the Best-Track observation file
fobs=open('BT_HISTORY_'+tcname,'r')
data_obs=numpy.loadtxt('BT_HISTORY_'+tcname,dtype='str')[:,:]

nbr_obs = 0
while fobs.readline():
    nbr_obs += 1
print('Il y a '+str(nbr_obs)+' positions Best-Track pour le cyclone '+tcname)

list_tc.append(tcname)


#### Loop on the number of analysed tropical cyclone
#
N=0  ### Initialize count for Invest or Nameless

for jname in range(nbr_tc):
    zradius=numpy.zeros(iradmax0, dtype=int)
    ztime=numpy.zeros(ech, dtype=object)
    ztime_basis=numpy.zeros(mb, dtype=object)

    iicentc=numpy.zeros(ech, dtype=int)
    ijcentc=numpy.zeros(ech, dtype=int)
    latcentc=numpy.zeros((mb,ech), dtype=float)
    loncentc=numpy.zeros((mb,ech), dtype=float)

    iicen=numpy.zeros(ech, dtype=object)
    ijcen=numpy.zeros(ech, dtype=object)

    zutcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zurcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zullcyl = numpy.zeros((ech,nphilS,iradmax0S), dtype=float)
    zvllcyl = numpy.zeros((ech,nphilS,iradmax0S), dtype=float)
    zuulcyl = numpy.zeros((ech,nphilS,iradmax0S), dtype=float)
    zvulcyl = numpy.zeros((ech,nphilS,iradmax0S), dtype=float)

    zrmwMEAN = numpy.zeros((ech,mb), dtype=float)

    zutcylMEANAZIM = numpy.zeros((ech,mb,iradmax0), dtype=float)
    zurcylMEANAZIM = numpy.zeros((ech,mb,iradmax0), dtype=float)

    zullcylmean = numpy.zeros((ech,mb,iradmax0S), dtype=float)
    zvllcylmean = numpy.zeros((ech,mb,iradmax0S), dtype=float)
    zuulcylmean = numpy.zeros((ech,mb,iradmax0S), dtype=float)
    zvulcylmean = numpy.zeros((ech,mb,iradmax0S), dtype=float)

    zddshear=numpy.zeros((ech,mb), dtype=float)
    zddsheardeg=numpy.zeros((ech,mb), dtype=float)
    zffshear=numpy.zeros((ech,mb), dtype=float)

    zwindsurfXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)

    zwindsurfXYQ90 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zwindsurfXYQ75 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zwindsurfXYQ50 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zwindsurfXYQ25 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)


    # Begin of the loop on member jb

    for jb in range(mb):
        print('On traite le membre MB='+str(list_mb[jb]).zfill(2))
        # Begin of the loop on time jc
        for jc in range(ech):
            print('On traite l\'échéance '+str(list_ech[jc]).zfill(2))
            #lecture  
            f = epygram.formats.resource('grid.arome-forecast.'+domaine+'+00'+str(list_ech[jc]).zfill(2)+':00-mb'+str(list_mb[jb]).zfill(3)+'.grib', 'r', 'GRIB')
            if(list_ech[jc]>0):
              fprec = epygram.formats.resource('grid.arome-forecast.'+domaine+'+00'+str(list_ech[jc-1]).zfill(2)+':00-mb'+str(list_mb[jb]).zfill(3)+'.grib', 'r', 'GRIB')  
            if (LGRIB2==True and list_mb[jb]<=16):
               mslp = f.readfield({'name': 'Pressure reduced to MSL'})
               ull=f.readfield({'name':'U component of wind', 'scaledValueOfFirstFixedSurface': 1500, 'level': 1500})
               vll=f.readfield({'name':'V component of wind', 'scaledValueOfFirstFixedSurface': 1500, 'level': 1500})
               uul=f.readfield({'name':'U component of wind', 'scaledValueOfFirstFixedSurface': 12000, 'level': 12000})
               vul=f.readfield({'name':'V component of wind', 'scaledValueOfFirstFixedSurface': 12000, 'level': 12000})
               u10=f.readfield({'name':'10 metre U wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
               v10=f.readfield({'name':'10 metre V wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
               mod10=numpy.sqrt(u10.data**2. + v10.data**2.)
            else:
               mslp = f.readfield({'indicatorOfParameter': 2})
            mslp.operation('/',100.)
            ztime_basis[jb]=mslp.validity.getbasis()
            ztime_term=mslp.validity.term()
            ztime[jc]=mslp.validity[0].get()
            print('Date traitee : '+str(ztime[jc]).zfill(2))
            #
            #
            # XX) Tracker initialization
            #
            if(list_ech[jc]%(6/echstep)==0):
              latobs, lonobs, vmaxobs, pminobs = read_BT(ficobs,tcname,ztime[jc])
            
            #
            #
            # XX) Finding TC's center
            #
            iicen[jc], ijcen[jc], lat, lon, latpos_1, lonpos_1, pmin = tracker_run(bassin,tcname,latobs,lonobs,speedguess,stepboxmove,xboxwind,mslp,u10,v10,iie,ije,zdeltax,zdeltay,list_ech,jc,jc-1)

            latcentc[jb,jc] = lat
            loncentc[jb,jc] = lon
            #
            #
            # XX) Calculating polar fields (u,v and magnitude)

            zutcyl[jc,:,:], zurcyl[jc,:,:] = polar_calcwind(mslp,u10.data,v10.data,iie,ije,iicen[jc],ijcen[jc],iradmax0,nphil,zdeltax,zdeltay,zdeltar,zdphi,Hemisphere)

            # XX) Azimuthally averaging

            zutcylMEANAZIM[jc,jb,:] = mean_azim(zutcyl[jc,:,:],iradmax0,nphil,flag=-999.)


            # XX) Compute the RMW

            zrmwMEAN[jc,jb] = compute_RMW(zutcylMEANAZIM[jc,jb,:],iradmax0,zdeltar)


            # XX) Computing wind shear

            zullcyl[jc,:,:] = polar_calc(mslp,ull.data,iie,ije,iicen[jc],ijcen[jc],iradmax0S,nphilS,zdeltax,zdeltay,zdeltarS,zdphiS)
            zvllcyl[jc,:,:] = polar_calc(mslp,vll.data,iie,ije,iicen[jc],ijcen[jc],iradmax0S,nphilS,zdeltax,zdeltay,zdeltarS,zdphiS)
            zuulcyl[jc,:,:] = polar_calc(mslp,uul.data,iie,ije,iicen[jc],ijcen[jc],iradmax0S,nphilS,zdeltax,zdeltay,zdeltarS,zdphiS)
            zvulcyl[jc,:,:] = polar_calc(mslp,vul.data,iie,ije,iicen[jc],ijcen[jc],iradmax0S,nphilS,zdeltax,zdeltay,zdeltarS,zdphiS)

            zullcylmean[jc,jb,:] = mean_azim(zullcyl[jc,:,:],iradmax0S,nphilS,flag=-999.)
            zvllcylmean[jc,jb,:] = mean_azim(zvllcyl[jc,:,:],iradmax0S,nphilS,flag=-999.)
            zuulcylmean[jc,jb,:] = mean_azim(zuulcyl[jc,:,:],iradmax0S,nphilS,flag=-999.)
            zvulcylmean[jc,jb,:] = mean_azim(zvulcyl[jc,:,:],iradmax0S,nphilS,flag=-999.)


            zffshear[jc,jb], zddsheardeg[jc,jb] = compute_windshear(zullcylmean[jc,jb,:],zvllcylmean[jc,jb,:],zuulcylmean[jc,jb,:],zvulcylmean[jc,jb,:],iradmax0S,zdeltarS)

            #
            ################################################################
            ### Compute Storm-relative wind (Vt, Ur) --> x,y
            ###############################################################
            #
            for jx in range(max(iicen[jc]-NX,1),min(iicen[jc]+NX+1,iie)):
                for jy in range(max(ijcen[jc]-NY,1),min(ijcen[jc]+NY+1,ije)):
                    zwindsurfXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = mod10.data[jy,jx] * 1.944      # Conversion en kts


zwindsurfXYQ90 = numpy.percentile(zwindsurfXY[:,:,:,:],90,axis=0)
zwindsurfXYQ75 = numpy.percentile(zwindsurfXY[:,:,:,:],75,axis=0)
zwindsurfXYQ50 = numpy.percentile(zwindsurfXY[:,:,:,:],50,axis=0)
zwindsurfXYQ25 = numpy.percentile(zwindsurfXY[:,:,:,:],25,axis=0)

# XX ) Mapping the storm-relative fields
for jc in range(ech):
    #
    #zoom=zoom-4.
    #zoom=zoom-4.
    centre_ll=mslp.geometry.ll2ij(numpy.mean(loncentc[:,jc]),numpy.mean(latcentc[:,jc]))
    centrectrl_ll=mslp.geometry.ll2ij(loncentc[0,jc],latcentc[0,jc])
    ijcentc[jc]=int(centre_ll[1])
    iicentc[jc]=int(centre_ll[0])

    #lattc, lontc = fix_TC_loc(ficEPS,mb,ztime[jc])
    locate_tc_wind_2D(u10,zwindsurfXY[0,jc,:,:],'CTRL',latcentc[0,jc],loncentc[0,jc],int(centrectrl_ll[0]),int(centrectrl_ll[1]),iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep,verifBT='true')
    locate_tc_wind_2D(u10,zwindsurfXYQ75[jc,:,:],'Q75',numpy.mean(latcentc[:,jc]),numpy.mean(loncentc[:,jc]),iicentc[jc],ijcentc[jc],iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep,verifBT='true')
    #
    im2 = Image.open('CTRL-WIND-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png')
    im3 = Image.open('Q75-WIND-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png')

    im2_size = im2.size
    im3_size = im3.size

    im_new = Image.new('RGB',(2*im2_size[0], im2_size[1]), (250,250,250))
    im_new.paste(im2,(0,0))
    im_new.paste(im3,(im2_size[0],0))


    #im_new = Image.new('RGB',(2*im1_size[0], 2*im1_size[1]), (250,250,250))
    #im_new.paste(im1,(0,0))
    #im_new.paste(im2,(im1_size[0],0))
    #im_new.paste(im3,(0,im1_size[0]))
    #im_new.paste(im4,(im1_size[0],im1_size[0]))

    im_new.save('VIEW-WIND-STORM-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png',"PNG")




# End of the loop on cyclones tcname      


