#!/usr/bin/env python
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
from plot_TC_features_tools import locate_tc_rain_2D

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
latn=float(data_nam[6])
lats=float(data_nam[7])
lono=float(data_nam[8])
lone=float(data_nam[9])
bassin=data_nam[10]
tcname=data_nam[11]

LGRIB2=True

zdeltax = 2500.
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
  iie=1395
  ije=747

if(domaine=='ncaled0025'):
  iie=521
  ije=491  


Alltrack=[]
Alltrackmed=[]
Alldate=[]
Allshear=[]
AnaPos=[]
AllTrace=[]
All64wind=[]
All34wind=[]
AllObs=[]
AllTimepos=[]

iicenG=572
ijcenG=252

speedguess = 50.
xradguess = 500. # en km
nphil=33  # number of azimuts
xboxwind=30. # 
xboxwindz=100.
iib=0
ijb=0

latinc=4
loninc=4

NX=80
NY=80

NXS=20
NYS=20

echdeb=0
echfin=75
echstep=3
stepboxmove = echstep


zdeltar = max(zdeltax,zdeltay)
zdphi = numpy.math.pi / 16.


#ficobs='BT_HISTORY_'+tcname+'_FAKE'
ficobs='BT_HISTORY_'+tcname
#ficEPS='PEindien-'+annee+mois+jour+heure+'-'+tcname

list_mb=[]
for m in range(0,16,1):
    list_mb.append(m) # Les niveaux lus dans les gribs Arome
mb=len(list_mb)

#list_levels=[500,2000,4000,6000,8000,10000,11000,12000,13000,14000,15000]
list_levels=[2000]
levels=len(list_levels)

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

    iicentc=numpy.zeros(ech, dtype=object)
    ijcentc=numpy.zeros(ech, dtype=object)
    #latobs=numpy.zeros(ech, dtype=object)
    #lonobs=numpy.zeros(ech, dtype=object)

    latcentc=numpy.zeros((mb,ech), dtype=float)
    loncentc=numpy.zeros((mb,ech), dtype=float)

    iicen=numpy.zeros(ech, dtype=object)
    ijcen=numpy.zeros(ech, dtype=object)

    zrain1XY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zrain2XY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zprecipXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)
    zprecipXY = numpy.zeros((mb,ech,2*NX+1,2*NY+1), dtype=float)

    zprecipXYQ90 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zprecipXYQ75 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zprecipXYQ50 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)
    zprecipXYQ25 = numpy.zeros((ech,2*NX+1,2*NY+1), dtype=float)

    X=numpy.zeros(2*NX+1, dtype=float)
    Y=numpy.zeros(2*NY+1, dtype=float)
    Xstar=numpy.zeros(2*NX+1, dtype=float)
    Ystar=numpy.zeros(2*NY+1, dtype=float)

    for jx in range(-NX,NX+1):
        X[jx+NX] = jx * zdeltax / 1000.
    for jy in range(-NY,NY+1):
        Y[jy+NY] = jy * zdeltay / 1000.


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
            #f = epygram.formats.resource('grid.arome-forecast.'+domaine+'+00'+str(list_ech[jc]).zfill(2)+':00.grib', 'r', 'GRIB')
            if (LGRIB2==True and list_mb[jb]<=17):
               mslp = f.readfield({'name': 'Pressure reduced to MSL'})
               u10=f.readfield({'name':'10 metre U wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
               v10=f.readfield({'name':'10 metre V wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
               mod10=numpy.sqrt(u10.data**2. + v10.data**2.)
               if(list_ech[jc]>0):
                 rain2 = f.readfield({'parameterCategory':1,'parameterNumber':65})
                 if(list_ech[jc]==echstep):
                    rain1 = fprec.readfield({'name': 'Pressure reduced to MSL'})
                    rain1.data[:,:] = 0.
                 else:
                    rain1 = fprec.readfield({'parameterCategory':1,'parameterNumber':65})
               else:
                   rain2 = f.readfield({'name': 'Pressure reduced to MSL'})
                   rain2.data[:,:] = 0.
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
            ################################################################
            ### Compute Storm-relative wind (Vt, Ur) --> x,y
            ###############################################################

            for jx in range(max(iicen[jc]-NX,1),min(iicen[jc]+NX+1,iie)):
                for jy in range(max(ijcen[jc]-NY,1),min(ijcen[jc]+NY+1,ije)):
                    if(list_ech[jc]>0):
                     zrain1XY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = rain1.data[jy,jx]
                     zrain2XY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = rain2.data[jy,jx]
                     zprecipXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zrain2XY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] - zrain1XY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]
                    else:
                     zrain2XY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = rain2.data[jy,jx]  
                     zprecipXY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX] = zrain2XY[jb,jc,jy-ijcen[jc]+NY,jx-iicen[jc]+NX]

#
zprecipXYQ90=numpy.percentile(zprecipXY[:,:,:,:],90,axis=0)
zprecipXYQ75=numpy.percentile(zprecipXY[:,:,:,:],75,axis=0)
zprecipXYQ50=numpy.percentile(zprecipXY[:,:,:,:],50,axis=0)
zprecipXYQ25=numpy.percentile(zprecipXY[:,:,:,:],25,axis=0)

# XX ) Mapping the storm-relative fields 
for jc in range(ech):
    #
    #zoom=zoom-4.
    #zoom=zoom-4.
    #centre_ll=mslp.geometry.ll2ij(numpy.mean(loncentc[:,jc]),numpy.mean(latcentc[:,jc]))
    centrectrl_ll=mslp.geometry.ll2ij(loncentc[0,jc],latcentc[0,jc])
    #ijcentc[jc]=int(centre_ll[1])
    #iicentc[jc]=int(centre_ll[0])
    #centre_ll=mslp.geometry.ll2ij(lonobs[jc],latobs[jc])

    #lattc, lontc = fix_TC_loc(ficEPS,mb,ztime[jc])
    locate_tc_rain_2D(u10,zprecipXY[0,jc,:,:],'CTRL',latcentc[0,jc],loncentc[0,jc],int(centrectrl_ll[0]),int(centrectrl_ll[1]),iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep)
    #locate_tc_rain_2D(u10,zprecipXY[0,jc,:,:],'CTRL',latobs[jc],lonobs[jc],int(centre_ll[0]),int(centre_ll[1]),iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep)
    #locate_tc_rain_2D(u10,zprecipXYQ75[jc,:,:],'Q75',numpy.mean(latcentc[:,jc]),numpy.mean(loncentc[:,jc]),iicentc[jc],ijcentc[jc],iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep)
    #locate_tc_rain_2D(u10,zprecipXYQ90[jc,:,:],'Q90',numpy.mean(latcentc[:,jc]),numpy.mean(loncentc[:,jc]),iicentc[jc],ijcentc[jc],iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep)
    #locate_tc_rain_2D(u10,zprecipXYQ90[jc,:,:],'Q90',latobs[jc],lonobs[jc],int(centre_ll[0]),int(centre_ll[1]),iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep)
    locate_tc_rain_2D(u10,zprecipXYQ90[jc,:,:],'Q90',latcentc[0,jc],loncentc[0,jc],int(centrectrl_ll[0]),int(centrectrl_ll[1]),iie,ije,NX,NY,zoom,ztime[jc],tcname,list_ech[jc],echstep)
    #
    im2 = Image.open('CTRL-RAIN-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png')
    im3 = Image.open('Q90-RAIN-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png')

    #im1_size = im1.size
    im2_size = im2.size
    im3_size = im3.size
    #im4_size = im4.size

    im_new = Image.new('RGB',(2*im2_size[0], im2_size[1]), (250,250,250))
    im_new.paste(im2,(0,0))
    im_new.paste(im3,(im2_size[0],0))

    #im_new = Image.new('RGB',(2*im1_size[0], 2*im1_size[1]), (250,250,250))
    #im_new.paste(im1,(0,0))
    #im_new.paste(im2,(im1_size[0],0))
    #im_new.paste(im3,(0,im1_size[0]))
    #im_new.paste(im4,(im1_size[0],im1_size[0]))

    im_new.save('VIEW-RR-STORM-CENTERED-'+tcname+'-'+str(list_ech[jc]).zfill(4)+'.png',"PNG")




# End of the loop on cyclones tcname      


