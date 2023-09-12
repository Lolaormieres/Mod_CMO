#!/usr/bin/env python
#!/bin/env python

from Magics.macro import *
import epygram
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

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
bassin=data_nam[6]
xp=data_nam[7]
tcname=data_nam[8]
echmax=int(data_nam[9])

LGRIB2=True

zdeltax = 2500.
zdeltay = 2500.

levelwind=1500


print("On se trouve sur le bassin=", bassin)

if(domaine=='GLOB01'):
  zdeltax = 10000.
  zdeltay = 10000.
  iie=3600
  ije=1801

if(domaine=='GLOB025'):
  zdeltax = 25000.
  zdeltay = 25000.
  iie=1440
  ije=721

if(domaine=='caraib0025'):
  iie=945
  ije=529

if(domaine=='INDIEN0025'):
  iie=1395
  ije=747

if(domaine=='ncaled0025'):
  iie=521
  ije=491  

if(domaine=='atllarge0025'):
  iie=1201
  ije=841

list_ech=[]
echstep=6
for t in range(0,echmax+echstep,echstep):
    list_ech.append(t) # Les niveaux lus dans les gribs Arome
ech=len(list_ech)
print('On a '+str(ech)+' echeances a traiter')


list_tc=[]


# Quelques parametres pour les calculs en geometrie cylindrique
#
speedguess = 50. # Vitesse d'advection moyenne du TC pour estimation position a l'instant t
xradguess = 200. # Radiale maximale en km
nphil=33  # Nombre d'azimuts (pour la projection cyclonique)
xboxwind=30. # Taille de la boite en kms pour la recherche du centre du TC (min de vent) 
zspeedmin=100. # Init pour le calcul du min de vent
iib=0
ijb=0

zdeltar = max(zdeltax,zdeltay)  # Resolution radiale (en prend la meme que la resolution de la grille)
zdphi = numpy.math.pi / 16.     # Resolution azimutale (ici pi/16.)

iradmax0 = numpy.int64((xradguess*1000. / zdeltar) + 2)  # Nombre de pts sur la radiale R


### For Visualisation
radmax     = xradguess


if(bassin=='CUSTOM-HN'):
   Hemisphere = 'N' 

if(bassin=='CUSTOM-HS'):
   Hemisphere = 'S' 
   
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
    ztime_basis=numpy.zeros(1,dtype=object)

    zutcylMEANAZIM = numpy.zeros((ech,iradmax0), dtype=float)
    zurcylMEANAZIM = numpy.zeros((ech,iradmax0), dtype=float)

    zucyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)    # Vent zonal U en (r,phi)
    zvcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)    # Vent meridien V en (r,phi)
    zutcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)   # Vent tangentiel en (r,phi)
    zurcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)   # Vent radial en (r,phi)

    #
    ## Write the data file for scores
    if os.path.exists(model+'-'+annee+mois+jour+heure+'-'+tcname):
       os.remove(model+'-'+annee+mois+jour+heure+'-'+tcname)
    fout=open(model+'-'+annee+mois+jour+heure+'-'+tcname,'a')
#
    # Begin of the loop on time jc  
    for jc in range(ech):
        zspeedmin=100.
        #lecture  
        f = epygram.formats.resource('grid.'+model+'-forecast.'+domaine+'+'+str(list_ech[jc]).zfill(4)+':00.grib', 'r', 'GRIB')
        if (LGRIB2==True):
           mslp = f.readfield({'name': 'Pressure reduced to MSL'})
        else:
           mslp = f.readfield({'indicatorOfParameter': 2})
        mslp.operation('/',100.)
        #
        ztime_basis[0]=mslp.validity.getbasis()
        ztime_term=mslp.validity.term()
        ztime[jc]=mslp.validity[0].get()
        for jobs in range(nbr_obs):
            if(str(data_obs[jobs,1]+' '+data_obs[jobs,2])==str(ztime[jc])):
               latobs=float(data_obs[jobs,3])
               lonobs=float(data_obs[jobs,4])
               vmaxobs=float(data_obs[jobs,5]) * 0.514  # conversion noeud --> m/s
               pminobs=float(data_obs[jobs,6])
               #latobs=float(data_obs[jobs,4])
               #lonobs=float(data_obs[jobs,5])
               #vmaxobs=float(data_obs[jobs,7]) * 0.514  # conversion noeud --> m/s
               #pminobs=float(data_obs[jobs,8])
        if(latobs==False):
           break
        ###############################################################
        # Ici on lit les parametres souhaites pour le calcul en
        # coordonnees polaires (r,phi)
        ###############################################################
        if (LGRIB2==True):
            u=f.readfield({'name':'U component of wind', 'scaledValueOfFirstFixedSurface': levelwind, 'level': levelwind})
            v=f.readfield({'name':'V component of wind', 'scaledValueOfFirstFixedSurface': levelwind, 'level': levelwind})
        else:
           u=f.readfield({'indicatorOfTypeOfLevel': 100, 'level': levelwind, 'table2Version': 1, 'editionNumber': 1, 'indicatorOfParameter': 33})
           v=f.readfield({'indicatorOfTypeOfLevel': 100, 'level': levelwind, 'table2Version': 1, 'editionNumber': 1, 'indicatorOfParameter': 34})
        
        ##########################################################################
        # Détermination de la boite de recherche des des conditions cycloniques 
        # (on considère un TC se déplaçant à une vitesse regiliere speedguess
        ##########################################################################
        if(jc-1>=0):
           data=numpy.loadtxt('pos.geo'+str(list_ech[jc-1]).zfill(2)+'-'+tcname)[:]
           latpos_1=data[0]
           lonpos_1=data[1]
        else:
           latguess=latobs
           longuess=lonobs
           if(bassin=='CUSTOM-HN'):
             print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % latguess)+' °N et '+str("% 4.2f" % (longuess*-1.))+' °W')
           if(bassin=='CUSTOM-HS'):
             print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % (latguess*-1.))+' °S et '+str("% 4.2f" % longuess)+' °E')
           latpos_1=latguess
           lonpos_1=longuess
        distguess = (echstep * speedguess) / 111.
        #
        latnpos = latpos_1 + distguess
        latspos = latpos_1 - distguess
        lonepos = lonpos_1 + distguess
        lonopos = lonpos_1 - distguess
        #
        niboxinf = max(int(mslp.geometry.ll2ij(lonopos,latspos)[0]),iib)
        niboxsup = min(int(mslp.geometry.ll2ij(lonepos,latspos)[0]),iie)
        njboxinf = max(int(mslp.geometry.ll2ij(lonopos,latspos)[1]),ijb)
        njboxsup = min(int(mslp.geometry.ll2ij(lonopos,latnpos)[1]),ije)
        #
        #
        if(niboxinf>=iib and niboxinf<=iie and njboxinf>=iib and njboxinf<=ije):
           #####################################################################
           # Détermination du centre cyclonique iicenpmin et ijcenpmin 
           # (minimum de pression réduite au niveau de la mer)
           #####################################################################
           mslpmin=1050
           for j in range(njboxinf,njboxsup):
               for i in range(niboxinf,niboxsup):
                   if mslp.data[j,i]<mslpmin:
                      mslpmin=mslp.data[j,i]
                      iicenpmin = i
                      ijcenpmin = j
#               
           centre_ll=mslp.geometry.ij2ll(iicenpmin,ijcenpmin)
#   
           #####################################################################
           # Détermination du centre cyclonique iicenvmin et ijcenvmin 
           # (minimum de vent dans l'oeil)
           #####################################################################
#
           if (LGRIB2==True):
              u10=f.readfield({'name':'10 metre U wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
              v10=f.readfield({'name':'10 metre V wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
           else:
              u10=f.readfield({'indicatorOfTypeOfLevel': 105, 'level': 10, 'table2Version': 1, 'editionNumber': 1, 'indicatorOfParameter': 33})
              v10=f.readfield({'indicatorOfTypeOfLevel': 105, 'level': 10, 'table2Version': 1, 'editionNumber': 1, 'indicatorOfParameter': 34})
#
           iavgwindi = int(round(xboxwind*1000./zdeltax))
           iavgwindj = int(round(xboxwind*1000./zdeltay))
#
           for jj in range(max(ijcenpmin-iavgwindj,1),min(ijcenpmin+iavgwindj,ije)):
               for ji in range(max(iicenpmin-iavgwindi,1),min(iicenpmin+iavgwindi,iie)):
                   zspeed=numpy.sqrt(u10.data[jj,ji]**2. + v10.data[jj,ji]**2.)
                   if zspeed<= zspeedmin:
                      zspeedmin=zspeed
                      iicenvmin = ji
                      ijcenvmin = jj
#  
#
           centre_ll=mslp.geometry.ij2ll(iicenvmin,ijcenvmin)
           if(bassin=='CUSTOM-HN'):
              print('A l\'échéance '+str(list_ech[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % centre_ll[1])+' °N et'+str("% 5.2f" % ((centre_ll[0]-360.)*-1)+' °W'))
           if(bassin=='CUSTOM-HS'):
              print('A l\'échéance '+str(list_ech[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % (centre_ll[1]*-1.))+' °S et'+str("% 5.2f" % centre_ll[0]+' °E'))   
#
           #####################################################################
           # Calcul des paramètres en géométrie cylindrique
           # (R=1 corresponds au centre du cyclone)
           #####################################################################
#           
           for jr in range(1,iradmax0):
               zradius[jr] = jr * zdeltar
# 
           zxi0=mslp.geometry.ij2xy(iicenvmin,ijcenvmin)[0]*100000. 
           zyj0=mslp.geometry.ij2xy(iicenvmin,ijcenvmin)[1]*100000. 
           zx00=mslp.geometry.ij2xy(1,ijcenvmin)[0]*100000. 
           zy00=mslp.geometry.ij2xy(iicenvmin,1)[1]*100000. 
#  
           for jr in range(iradmax0):
               for jphi in range(nphil):
                   zphi = (jphi) * zdphi
                   zxk  = zradius[jr] * numpy.cos(zphi) + zxi0
                   zyk  = zradius[jr] * numpy.sin(zphi) + zyj0
                   iix  = int(((zxk-zx00) / zdeltax) + 1 )
                   iiy  = int(((zyk-zy00) / zdeltay) + 1)
#
                   if (iix>iib) and ((iix+1)<iie) and (iiy>ijb) and ((iiy+1)<ije): 
                      zxk = (zxk-mslp.geometry.ij2xy(iix,iiy)[0]*100000.) / zdeltax 
                      zyk = (zyk-mslp.geometry.ij2xy(iix,iiy)[1]*100000.) / zdeltay 
                      zucyl[jc,jphi,jr]=u.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                        u.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                        u.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                        u.data[iiy+1,iix+1]*zxk*zyk
                      zvcyl[jc,jphi,jr]=v.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                        v.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                        v.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                        v.data[iiy+1,iix+1]*zxk*zyk
#
                      if Hemisphere=='S':
                         zucyl[jc,jphi,jr] =  -zucyl[jc,jphi,jr]
                         zvcyl[jc,jphi,jr] =  -zvcyl[jc,jphi,jr]
                      # Calcul du vent tangentiel Vt   
                      zutcyl[jc,jphi,jr] = zvcyl[jc,jphi,jr]*numpy.cos(zphi) \
                                           - zucyl[jc,jphi,jr]*numpy.sin(zphi)
                      zurcyl[jc,jphi,jr] = zucyl[jc,jphi,jr]*numpy.cos(zphi) \
                                                      + zvcyl[jc,jphi,jr]*numpy.sin(zphi)                                           
                      if zutcyl[jc,jphi,jr]<0: 
                         zutcyl[jc,jphi,jr] = 0.
                   #
                   else:
                      zutcyl[jc,jphi,jr]=-999. 
                      zurcyl[jc,jphi,jr]=-999.
               icount=0
            #
            #
               icount=0
               for jphi in range(0,32):
                   if(zutcyl[jc,jphi,jr]!=-999.):
                      zutcylMEANAZIM[jc,jr] = zutcylMEANAZIM[jc,jr] + zutcyl[jc,jphi,jr]
                      zurcylMEANAZIM[jc,jr] = zurcylMEANAZIM[jc,jr] + zurcyl[jc,jphi,jr]
                      icount = icount + 1
                   else:
                      zutcylMEANAZIM[jc,jr]=-999.
                      zurcylMEANAZIM[jc,jr]=-999.
               if icount == len(range(0,32)):
                  zutcylMEANAZIM[jc,jr] = zutcylMEANAZIM[jc,jr] / icount
                  zurcylMEANAZIM[jc,jr] = zurcylMEANAZIM[jc,jr] / icount
               else:
                  zutcylMEANAZIM[jc,jr] = -999.
                  zurcylMEANAZIM[jc,jr] = -999.
            #
            #
#       
           # On recupere les lat/lon du TC relocalise a l'instant t
           lat=centre_ll[1]
           lon=centre_ll[0]-360.
           vmax=numpy.max(zutcyl[jc,:,:])
           pmin=mslpmin

           if os.path.exists('pos.geo'+str(list_ech[jc]).zfill(2)+'-'+tcname):
              os.remove('pos.geo'+str(list_ech[jc]).zfill(2)+'-'+tcname)
           f=open('pos.geo'+str(list_ech[jc]).zfill(2)+'-'+tcname,'a')
           f.write("#GEO\n")
           f.write("#FORMAT XY_VECTOR\n")
           f.write("# lat lon\n")
           f.write("#DATA\n")
           f.write(str("%4.2f" % lat)+' '+str("%4.2f" % lon)+"\n")
           f.close()
           #
           # Ecriture des positions du TC a chaque echeance
           fout.write(tcname+' '+str(ztime[jc])+' '+str("%4.2f" % lat)+' '+str("%4.2f" % lon)+' '+str("%4.1f" % vmax)+' '+str("%4.1f" % pmin)+"\r\n")
##
#### Trace des champs moyennes sur l'azimut en fonction du temps (Hovmoller)
#
    zradius=zradius/1000. # En kms

    cmap1 = mpl.colors.ListedColormap([(0.9,0.9,0.9),(0.7,0.7,0.7),(0.5,0.5,0.5),(0.3,0.3,0.3),(0.1,0.1,0.1)])
    cmap1.set_over('black')
    cmap1.set_under('white')
#
    fig, a=plt.subplots(nrows=1, ncols=1, figsize=(25,15))
#
    vt= a.contourf(zradius[:],ztime[:],zutcylMEANAZIM[:,:],levels=[0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80], cmap='jet')
    #ur=a.contour(zradius[:],ztime[:],zurcylMEANAZIM[:,:], levels=[-6,-4,-2], colors='black',linewidths=3,linestyles='solid')
    #a.clabel(ur,fmt='%3i',colors='black',fontsize=20)
#
    a.contourf(zradius[:],ztime[:],zutcylMEANAZIM[:,:],1, colors='None',hatches=['O'],levels=[-999,0])
#
    #cbar=fig.colorbar(ur.lines,ticks=[0.,5.,10.,15.,20.,25.])
    #cbar.set_clim(0,25)
    #cbar.draw_all()
    cbar=fig.colorbar(vt, ax=a)
    a.set_title(tcname+':Hovmöller Vt+Ur @ '+str(levelwind)+' m - Time : '+str(ztime[jc]), loc='right', fontsize=32, fontstyle='normal')
    a.set_xlabel('Distance depuis le centre (km)', fontsize=48)
    a.set_ylabel('Heure (UTC)', fontsize=48)
    a.set_xlim(0,radmax)
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    cbar.ax.tick_params(labelsize=42)
    cbar.set_label('(m/s)',fontsize=66)

    fig.savefig(xp+'_Hovmo-'+tcname+'-MeanAZIM-'+annee+mois+jour+heure+'testGLOB01.png',dpi=200)   

fout.close()
