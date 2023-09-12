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
from compute_TC_features_tools import polar_calc, polar_calcwind, mean_azim, read_BT, fix_TC_loc, write_model_data, tracker_run, compute_RMW, compute_Vtmax_azim, compute_windshear, compute_qv, compute_THETAE
from plot_TC_features_tools import locate_tc_wind_2D, locate_tc_rain_2D, locate_tc_isp_2D, locate_tc_innerproc_2D

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

if(domaine=='GLOB01'):
  Hemisphere = 'N'  
  zdeltax = 10000.
  zdeltay = 10000.
  iie=3600
  ije=1801
  print('*****************************laaaaaaaaaaaaaaaaaaaaaaaa***********************************')

if(domaine=='GLOB025'):
  Hemisphere = 'N'  
  zdeltax = 25000.
  zdeltay = 25000.
  iie=1440
  ije=721

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

if(domaine=='ncaled0025'):
  iie=521
  ije=491  



speedguess = 50. # Vitesse d'advection de la boîte pour la recherche du centre cyclonique
xradguess = 800. # Distance depuis le centre du cyclone en coordonnées polaires en km
zdeltar = max(zdeltax,zdeltay)    # Résolution radiale en polaire (même que la grille de sortie)
zdphi = numpy.math.pi / 16. # Résolution azimutale
nphil = int((2. * numpy.math.pi / zdphi)) + 1 # Nbre d'azimuts en fonction de la résolution
xboxwind=30. # Boite autour de la Pmin pour trouver le min de vent comme centre (attention désactivé pour le moment dans tracker_run)
iib=0
ijb=0
seuildbz=30.

nivzmin=500.
nivzmax=2000.
nivzstep=500.
nstepsel=numpy.int64((nivzmax-nivzmin)/nivzstep + 1)

KMX=400.   # Distance souhaitée autour du centre cyclonique dans le repère cartésien (en km)
KMY=400.

NX=int(KMX/(zdeltax/1000.))    
NY=int(KMY/(zdeltay/1000.))

#   Fixer les échéances et pas de temps ici
echdeb=0
echfin=echmax+6
echstep=6
stepboxmove = echstep



ficobs='BT_HISTORY_'+tcname

list_mb=[]
for m in range(0,1,1):
    list_mb.append(m) # Les niveaux lus dans les gribs Arome
mb=len(list_mb)
print('mb',mb)

#list_levels=[500,2000,4000,6000,8000,10000,11000,12000,13000,14000,15000]
#list_levels=[500]
#levels=len(list_levels)

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
nphilS=73

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
    list_nivz=numpy.zeros(nstepsel, dtype=object)

    iicentc=numpy.zeros(ech, dtype=int)
    ijcentc=numpy.zeros(ech, dtype=int)
    iicentcfix=numpy.zeros(ech, dtype=int)
    ijcentcfix=numpy.zeros(ech, dtype=int)
    latcentc=numpy.zeros((mb,ech), dtype=float)
    loncentc=numpy.zeros((mb,ech), dtype=float)

    iicen=numpy.zeros(ech, dtype=object)
    ijcen=numpy.zeros(ech, dtype=object)

    zutcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zurcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zuurcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zvurcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)

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
    zurXY=numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zvtXY=numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zuflow=numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zvflow=numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zuurXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zvurXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zuvtXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zvvtXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zmodvtXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zwindsurfXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zdbz12kmXY  = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zmucapeXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)

    zwindsurfXYQ90 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zurXYQ10 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zwindsurfXYQ75 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zwindsurfXYQ50 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zwindsurfXYQ25 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zuurXYQ10 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zvurXYQ10 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zthaeXY=numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)

    zmaskcb1XY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zmaskcb1XYP = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)

    uradtempo = numpy.zeros((ije,iie), dtype=float)
    vradtempo = numpy.zeros((ije,iie), dtype=float)


    # Begin of the loop on member jb

    for jb in range(mb):
        print('On traite le membre MB='+str(list_mb[jb]).zfill(2))
        # Begin of the loop on time jc
        for jc in range(ech):
            print('On traite l\'échéance '+str(list_ech[jc]).zfill(2))
            #lecture  
            f = epygram.formats.resource('grid.arpege-forecast.'+domaine+'+00'+str(list_ech[jc]).zfill(2)+':00.grib', 'r', 'GRIB')
            #f = epygram.formats.resource('grid.arpege-forecast.'+domaine+'+00'+str(list_ech[jc]).zfill(2)+':00-mb'+str(list_mb[jb]).zfill(3)+'.grib', 'r', 'GRIB')
            if(list_ech[jc]>0):
              #fprec = epygram.formats.resource('grid.arpege-forecast.'+domaine+'+00'+str(list_ech[jc-1]).zfill(2)+':00-mb'+str(list_mb[jb]).zfill(3)+'.grib', 'r', 'GRIB')  
              fprec = epygram.formats.resource('grid.arpege-forecast.'+domaine+'+00'+str(list_ech[jc-1]).zfill(2)+':00.grib', 'r', 'GRIB')  
            if (LGRIB2==True):
               mslp = f.readfield({'name': 'Pressure reduced to MSL'})
               #ull=f.readfield({'name':'U component of wind', 'scaledValueOfFirstFixedSurface': 1500, 'level': 1500})
               #vll=f.readfield({'name':'V component of wind', 'scaledValueOfFirstFixedSurface': 1500, 'level': 1500})
               #uul=f.readfield({'name':'U component of wind', 'scaledValueOfFirstFixedSurface': 12000, 'level': 12000})
               #vul=f.readfield({'name':'V component of wind', 'scaledValueOfFirstFixedSurface': 12000, 'level': 12000})
               uul=f.readfield({'name':'U component of wind', 'scaledValueOfFirstFixedSurface' : 200.*100., 'level': 200.})
               vul=f.readfield({'name':'V component of wind', 'scaledValueOfFirstFixedSurface' : 200.*100., 'level': 200.})
               u10=f.readfield({'name':'10 metre U wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
               v10=f.readfield({'name':'10 metre V wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})

               mod10=numpy.sqrt(u10.data**2. + v10.data**2.)
               #dbz12km = f.readfield({'parameterNumber': 4, 'scaledValueOfFirstFixedSurface': 12000, 'level': 12000})

               #p500=f.readfield({'parameterNumber':0,'parameterCategory':3, 'scaledValueOfFirstFixedSurface': 500, 'level': 500})
               #t500=f.readfield({'parameterNumber':0,'parameterCategory':0, 'scaledValueOfFirstFixedSurface': 500, 'level': 500})
               #hu500=f.readfield({'parameterNumber':1,'parameterCategory':1, 'scaledValueOfFirstFixedSurface': 500, 'level': 500})

               #p5km=f.readfield({'parameterNumber':0,'parameterCategory':3, 'scaledValueOfFirstFixedSurface': +nivzmax, 'level': +nivzmax})
               #t5km=f.readfield({'parameterNumber':0,'parameterCategory':0, 'scaledValueOfFirstFixedSurface': +nivzmax, 'level': +nivzmax})
               #hu5km=f.readfield({'parameterNumber':1,'parameterCategory':1, 'scaledValueOfFirstFixedSurface': +nivzmax, 'level': +nivzmax})

               #zqv500 = compute_qv(p500,t500,hu500,iie,ije)
               #zqv5km = compute_qv(p5km,t5km,hu5km,iie,ije)
               #zqv500 = zqv500 * 1000.     # --> g/kg
               #zqv5km = zqv5km * 1000.     # --> g/kg

               #p500.operation('/',100.)
               #p5km.operation('/',100.)

               #thetae500 = compute_THETAE(p500,t500,zqv500,iie,ije)
               #thetae5km = compute_THETAE(p5km,t5km,zqv5km,iie,ije)

               urad=0.
               vrad=0.
               thae=0.
               ncount=0
               for z in range(nstepsel):
                   list_nivz[z]= nivzmin + z*nivzstep # Les niveaux lus dans les gribs Arome
                   print('On lit le niveau z='+str(list_nivz[z]).zfill(3)+' m')
                   urad=f.readfield({'name':'U component of wind', 'scaledValueOfFirstFixedSurface': +list_nivz[z], 'level': +list_nivz[z]})
                   vrad=f.readfield({'name':'V component of wind', 'scaledValueOfFirstFixedSurface': +list_nivz[z], 'level': +list_nivz[z]})

                   uradtempo=uradtempo+urad.data
                   vradtempo=vradtempo+vrad.data
                   ncount=ncount+1
               print('Les basses couches sont moyennées sur '+str(ncount)+' niveaux')
               if(ncount!=0):
                  uradtempo=uradtempo/ncount
                  vradtempo=vradtempo/ncount
                  urad.data=uradtempo
                  vrad.data=vradtempo
               #print("Valeur max de thetae=",numpy.max(thetae))

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

            #zutcyl[jc,:,:], zurcyl[jc,:,:] = polar_calcwind(mslp,u10.data,v10.data,iie,ije,iicen[jc],ijcen[jc],iradmax0,nphil,zdeltax,zdeltay,zdeltar,zdphi,Hemisphere)

            # XX) Azimuthally averaging

            #zutcylMEANAZIM[jc,jb,:] = mean_azim(zutcyl[jc,:,:],iradmax0,nphil,flag=-999.)


            # XX) Compute the RMW

            #zrmwMEAN[jc,jb] = compute_RMW(zutcylMEANAZIM[jc,jb,:],iradmax0,zdeltar)


            # XX) Computing wind shear

            #zullcyl[jc,:,:] = polar_calc(mslp,ull.data,iie,ije,iicen[jc],ijcen[jc],iradmax0S,nphilS,zdeltax,zdeltay,zdeltarS,zdphiS)
            #zvllcyl[jc,:,:] = polar_calc(mslp,vll.data,iie,ije,iicen[jc],ijcen[jc],iradmax0S,nphilS,zdeltax,zdeltay,zdeltarS,zdphiS)
            #zuulcyl[jc,:,:] = polar_calc(mslp,uul.data,iie,ije,iicen[jc],ijcen[jc],iradmax0S,nphilS,zdeltax,zdeltay,zdeltarS,zdphiS)
            #zvulcyl[jc,:,:] = polar_calc(mslp,vul.data,iie,ije,iicen[jc],ijcen[jc],iradmax0S,nphilS,zdeltax,zdeltay,zdeltarS,zdphiS)

            #zullcylmean[jc,jb,:] = mean_azim(zullcyl[jc,:,:],iradmax0S,nphilS,flag=-999.)
            #zvllcylmean[jc,jb,:] = mean_azim(zvllcyl[jc,:,:],iradmax0S,nphilS,flag=-999.)
            #zuulcylmean[jc,jb,:] = mean_azim(zuulcyl[jc,:,:],iradmax0S,nphilS,flag=-999.)
            #zvulcylmean[jc,jb,:] = mean_azim(zvulcyl[jc,:,:],iradmax0S,nphilS,flag=-999.)


            #zffshear[jc,jb], zddsheardeg[jc,jb] = compute_windshear(zullcylmean[jc,jb,:],zvllcylmean[jc,jb,:],zuulcylmean[jc,jb,:],zvulcylmean[jc,jb,:],iradmax0S,zdeltarS)

            #print("Max de zddsheardeg=",numpy.max(zddsheardeg[jc,jb]))
            #zddshear[jc,jb] = (zddsheardeg[jc,jb] + 0.) * numpy.math.pi / 180. #

            #print("Max de zddshear=",numpy.max(zddshear[jc,jb]))




            #
            #
            ################################################################
            ### Compute Storm-relative wind (Vt, Ur) --> x,y
            ###############################################################
            #
            for jx in range(max(iicen[jc]-NX,1),min(iicen[jc]+NX+1,iie)):
                for jy in range(max(ijcen[jc]-NY,1),min(ijcen[jc]+NY+1,ije)):
                    zwindsurfXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = mod10.data[jy,jx] * 1.944      # Conversion en kts
                    #zdbz12kmXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = dbz12km.data[jy,jx]
                    #zthaeXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]= thetae5km[jy,jx] - thetae500[jy,jx]
                    zuflow[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]= urad.data[jy,jx]
                    zvflow[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]= vrad.data[jy,jx]
                    #if(zdbz12kmXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] > seuildbz):
                    #  zmaskcb1XY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = 1.  
                    #
                    zxk=mslp.geometry.ij2xy(jx,ijcen[jc])[0]*100000. - mslp.geometry.ij2xy(iicen[jc],ijcen[jc])[0]*100000.
                    zyk=mslp.geometry.ij2xy(iicen[jc],jy)[1]*100000. - mslp.geometry.ij2xy(iicen[jc],ijcen[jc])[1]*100000.
                    zr2=zxk*zxk + zyk*zyk
                    zr=numpy.sqrt(zr2)
                    zphi = numpy.math.atan2(zyk,zxk)
                    zurXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zuflow[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]*numpy.cos(zphi)+zvflow[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]*numpy.sin(zphi)
                    zvtXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zvflow[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]*numpy.cos(zphi)-zuflow[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]*numpy.sin(zphi)
                    if(zurXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]>=0.):
                      zurXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]=0.
                    zuurXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zurXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] * numpy.cos(zphi)
                    zvurXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zurXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] * numpy.sin(zphi)

                    zuvtXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = - zvtXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] * numpy.sin(zphi)
                    zvvtXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zvtXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] * numpy.cos(zphi)

                    zmodvtXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = numpy.sqrt(zuvtXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]**2. + zvvtXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]**2.)



#zurXYQ10 = numpy.percentile(zurXY[:,:,:,:],10,axis=0)

#for jc in range(ech):
#    for jx in range(max(iicen[jc]-NX,1),min(iicen[jc]+NX+1,iie)):
#        for jy in range(max(ijcen[jc]-NY,1),min(ijcen[jc]+NY+1,ije)):
#            zxk=mslp.geometry.ij2xy(jx,ijcen[jc])[0]*100000. - mslp.geometry.ij2xy(iicen[jc],ijcen[jc])[0]*100000.
#            zyk=mslp.geometry.ij2xy(iicen[jc],jy)[1]*100000. - mslp.geometry.ij2xy(iicen[jc],ijcen[jc])[1]*100000.
#            zr2=zxk*zxk + zyk*zyk
#            zr=numpy.sqrt(zr2)
#            zphi = numpy.math.atan2(zyk,zxk)
#            if(zurXYQ10[jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]>=0.):
#              zurXYQ10[jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]=0.
#            zuurXYQ10[jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zurXYQ10[jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] * numpy.cos(zphi)
#            zvurXYQ10[jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zurXYQ10[jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] * numpy.sin(zphi)

zwindsurfXYQ90 = numpy.percentile(zwindsurfXY[:,:,:,:],90,axis=0)
zwindsurfXYQ75 = numpy.percentile(zwindsurfXY[:,:,:,:],75,axis=0)
zwindsurfXYQ50 = numpy.percentile(zwindsurfXY[:,:,:,:],50,axis=0)
zwindsurfXYQ25 = numpy.percentile(zwindsurfXY[:,:,:,:],25,axis=0)

#zuurXYQ90 = numpy.percentile(zuurXY[:,:,:,:],90,axis=0)
#zvurXYQ90 = numpy.percentile(zvurXY[:,:,:,:],90,axis=0)
#zuurXYQ10 = numpy.percentile(zuurXY[:,:,:,:],10,axis=0)
#zvurXYQ10 = numpy.percentile(zvurXY[:,:,:,:],10,axis=0)

#zthaeXYQ10=numpy.percentile(zthaeXY[:,:,:,:],10,axis=0)
#zthaeXYQ90=numpy.percentile(zthaeXY[:,:,:,:],90,axis=0)


#zmaskcb1XYP[:,:,:] = (numpy.sum(zmaskcb1XY[:,:,:,:],axis=0) / mb ) * 100.


# XX ) Mapping the storm-relative fields
for jc in range(ech):
    #
    #zmaskcb1XY[0,jc,:,:] = zmaskcb1XY[0,jc,:,:] * 99.
    #print("Max de zmaskcb1XYP=",numpy.max(zmaskcb1XYP[jc,:,:]))
    #zoom=zoom-4.
    #zoom=zoom-4.
    centre_ll=mslp.geometry.ll2ij(numpy.mean(loncentc[:,jc]),numpy.mean(latcentc[:,jc]))
    centrectrl_ll=mslp.geometry.ll2ij(loncentc[0,jc],latcentc[0,jc])
    ijcentc[jc]=int(centre_ll[1])
    iicentc[jc]=int(centre_ll[0])
    #
    locate_tc_wind_2D(mslp,zmodvtXY[0,jc,:,:],u10,v10,zuurXY[0,jc,:,:],zvurXY[0,jc,:,:],'OPER',latcentc[0,jc],loncentc[0,jc],int(centrectrl_ll[0]),int(centrectrl_ll[1]),iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep,verifBT='true')
    #locate_tc_wind_2D(u10,zwindsurfXYQ75[jc,:,:],'Q75',numpy.mean(latcentc[:,jc]),numpy.mean(loncentc[:,jc]),iicentc[jc],ijcentc[jc],iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep,verifBT='true')
    #
    #im2 = Image.open('OPER-WIND-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png')
    #im3 = Image.open('Q90-CB-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png')

    #im2_size = im2.size
    #im3_size = im3.size

    #im_new = Image.new('RGB',(2*im2_size[0], im2_size[1]), (250,250,250))
    #im_new.paste(im2,(0,0))
    #im_new.paste(im3,(im2_size[0],0))


    #im_new = Image.new('RGB',(2*im1_size[0], 2*im1_size[1]), (250,250,250))
    #im_new.paste(im1,(0,0))
    #im_new.paste(im2,(im1_size[0],0))
    #im_new.paste(im3,(0,im1_size[0]))
    #im_new.paste(im4,(im1_size[0],im1_size[0]))

    #im_new.save('VIEW-CB-STORM-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png',"PNG")




# End of the loop on cyclones tcname      


