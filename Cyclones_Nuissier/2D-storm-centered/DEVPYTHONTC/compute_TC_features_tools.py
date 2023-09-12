
def read_BT(ficBT,tcname,time):

    import numpy
    import os
    
    ### Read the Best-Track observation file
    fobs=open(ficBT,'r')
    data_obs=numpy.loadtxt(ficBT,dtype='str')[:,:]

    nbr_obs = 0
    while fobs.readline():
          nbr_obs += 1
    print('Il y a '+str(nbr_obs)+' positions Best-Track pour le cyclone '+tcname)

    dateOK=0
    for jobs in range(nbr_obs):
        if(str(data_obs[jobs,1]+' '+data_obs[jobs,2])==str(time)):
           latobs=float(data_obs[jobs,3])
           lonobs=float(data_obs[jobs,4])
           vmaxobs=float(data_obs[jobs,5]) * 0.514  # conversion noeud --> m/s
           pminobs=float(data_obs[jobs,6])
           #latobs=float(data_obs[jobs,4])
           #lonobs=float(data_obs[jobs,5])
           #vmaxobs=float(data_obs[jobs,7]) * 0.514  # conversion noeud --> m/s
           #pminobs=float(data_obs[jobs,8])
           #rmwobs=float(data_obs[jobs,19])  * 1.61  # conversion miles --> km
           dateOK=dateOK+1
    if(dateOK==0):
       print("Attention pas de date correspondante dans le fichier BT !!") 
       os.abort() 

    return latobs, lonobs, vmaxobs, pminobs

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def fix_TC_loc(ficEPS,mb,time):

    import numpy
    import os

    lateps = numpy.zeros(mb, dtype=float)
    loneps = numpy.zeros(mb, dtype=float)

    ### Read the EPS data file
    feps=open(ficEPS,'r')
    data_eps=numpy.loadtxt(ficEPS,dtype='str')[:,:]

    nbr_eps = 0
    while feps.readline():
          nbr_eps += 1

    for jb in range(mb):
        for jdata in range(nbr_eps):
            if((str(data_eps[jdata,2]+' '+data_eps[jdata,3])==str(time)) and (float(data_eps[jdata,4])==jb+1)):
              lateps[jb]=float(data_eps[jdata,5])
              loneps[jb]=float(data_eps[jdata,6])

    latloc=numpy.mean(lateps)
    lonloc=numpy.mean(loneps)

    return latloc, lonloc

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def write_model_data(fout,tcname,ztime,lat,lon):

    fout.write(tcname+' '+str(ztime)+' '+str("%4.2f" % lat)+' '+str("%4.2f" % lon)+"\r\n")

    return

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def tracker_run(bassin,tcname,latobs,lonobs,speedguess,stepboxmove,xboxwind,mslp,u10,v10,iie,ije,zdeltax,zdeltay,status,jc,jc_1):
    ## *TRACKER_RUN* - 
    ## 
    ## PURPOSE
    ## -------
    # This function does TC's tracking over the geographical domain using minimum central pressure (Pmin) and local wind minimum (Vmin)
    # status=list of time, jc stands for a given t time whearas jc_1 for a given t-1 time


    import numpy
    import os

    iib=0
    ijb=0


    if(jc_1>=0):
      data=numpy.loadtxt('pos.geo'+str(status[jc_1]).zfill(2)+'-'+tcname)[:]
      lat_1=data[0]
      lon_1=data[1]
    else:
      latguess=latobs
      longuess=lonobs
      lat_1=latguess
      lon_1=longuess
      if(bassin=='CUSTOM-HN'):
        print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % latguess)+' °N et '+str("% 4.2f" % (longuess*-1.))+' °W')
      if(bassin=='CUSTOM-HS'):
        print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % (latguess*-1.))+' °S et '+str("% 4.2f" % longuess)+' °E')
    distguess = (stepboxmove * speedguess) / 111.
    #
    latnpos = lat_1 + distguess
    latspos = lat_1 - distguess
    lonepos = lon_1 + distguess
    lonopos = lon_1 - distguess
    #
    niboxinf = max(int(mslp.geometry.ll2ij(lonopos,latspos)[0]),iib)
    niboxsup = min(int(mslp.geometry.ll2ij(lonepos,latspos)[0]),iie)
    njboxinf = max(int(mslp.geometry.ll2ij(lonopos,latspos)[1]),ijb)
    njboxsup = min(int(mslp.geometry.ll2ij(lonopos,latnpos)[1]),ije)
    #
    if(niboxinf>=iib and niboxinf<=iie and njboxinf>=iib and njboxinf<=ije):
      #####################################################################
      # Looking for the minimum central MSLP
      # (as a firts guess of the TC's center)
      #####################################################################
      mslpmin=1050
      for j in range(njboxinf,njboxsup):
          for i in range(niboxinf,niboxsup):
              if mslp.data[j,i]<mslpmin:
                 mslpmin=mslp.data[j,i]
                 iicen = i
                 ijcen = j
#
      #####################################################################
      # Looking for the minimum of local horizontal wind near Pmin
      # (as the fixed TC's center)
      #####################################################################
#
#
      #iavgwindi = int(round(xboxwind*1000./zdeltax))
      #iavgwindj = int(round(xboxwind*1000./zdeltay))
#
      #zspeedmin=100.
      #for jj in range(max(ijcen-iavgwindj,1),min(ijcen+iavgwindj,ije)):
      #    for ji in range(max(iicen-iavgwindi,1),min(iicen+iavgwindi,iie)):
      #        zspeed=numpy.sqrt(u10.data[jj,ji]**2. + v10.data[jj,ji]**2.)
      #        if zspeed < zspeedmin:
      #           zspeedmin=zspeed
      #           iicen = ji
      #           ijcen = jj

      #for jj in range(ijcen-iavgwindj,ijcen+iavgwindj):
      #    for ji in range(iicen-iavgwindi,iicen+iavgwindi):
      #        if mod10.data[jj,ji]<= zspeedmin:
      #           zspeedmin=mod10.data[jj,ji]
      #           iicen = ji
      #           ijcen = jj
#
      #vortmax=0.
      #for jj in range(max(ijcen-iavgwindj,1),min(ijcen+iavgwindj,ije)):
      #    for ji in range(max(iicen-iavgwindi,1),min(iicen+iavgwindi,iie)):
      #        if vortll.data[jj,ii]>vortmax:
      #           vortmax=vortll.data[jj,ii]
      #           iicen = ji
      #           ijcen = jj

      centre_ll=mslp.geometry.ij2ll(iicen,ijcen)
      lat=centre_ll[1]
      lon=centre_ll[0]-360.
      if(bassin=='CUSTOM-HN'):
        print('A l\'échéance '+str(status[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % centre_ll[1])+' °N et'+str("% 5.2f" % ((centre_ll[0]-360.)*-1)+' °W'))
      if(bassin=='CUSTOM-HS'):
        print('A l\'échéance '+str(status[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % (centre_ll[1]*-1.))+' °S et'+str("% 5.2f" % centre_ll[0]+' °E'))

      #
      if os.path.exists('pos.geo'+str(status[jc]).zfill(2)+'-'+tcname):
         os.remove('pos.geo'+str(status[jc]).zfill(2)+'-'+tcname)
      f=open('pos.geo'+str(status[jc]).zfill(2)+'-'+tcname,'a')
      f.write("#GEO\n")
      f.write("#FORMAT XY_VECTOR\n")
      f.write("# lat lon height date time u v\n")
      f.write("#DATA\n")
      f.write(str("%4.2f" % lat)+' '+str("%4.2f" % lon)+' 0 20030617 1200 '+str("%4.2f" % mslpmin)+' '+str("%4.2f" % mslpmin)+"\n")
      f.close()


    return iicen, ijcen, lat, lon, lat_1, lon_1, mslpmin

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def polar_calc(dom2readgeom,zvarcart,iie,ije,icen,jcen,iradmax0,nphil,zdeltax,zdeltay,zdeltar,zdphi):

    ## *POLAR_CALC* - interpolation onto cylindrical grid
    ## 
    ## PURPOSE
    ## -------
    # This subroutine interpolates an input field onto a cylindrical grid centred on (icen,jcen) at a given level (R=1 correspond to the TC's center)

    # INITIALIZATION
    #
    import numpy
    import epygram

    zradius=numpy.zeros(iradmax0, dtype=int)
    zvarcyl = numpy.zeros((nphil,iradmax0), dtype=float)

    iib=0
    ijb=0


    for jr in range(1,iradmax0):
        zradius[jr] = jr * zdeltar
    # 
    zxi0=dom2readgeom.geometry.ij2xy(icen,jcen)[0]*100000.
    zyj0=dom2readgeom.geometry.ij2xy(icen,jcen)[1]*100000.
    zx00=dom2readgeom.geometry.ij2xy(1,jcen)[0]*100000.
    zy00=dom2readgeom.geometry.ij2xy(icen,1)[1]*100000.

    for jr in range(iradmax0):
        for jphi in range(nphil):
            zphi = (jphi) * zdphi
            zxk  = zradius[jr] * numpy.cos(zphi) + zxi0
            zyk  = zradius[jr] * numpy.sin(zphi) + zyj0
            iix  = int(((zxk-zx00) / zdeltax) + 1 )
            iiy  = int(((zyk-zy00) / zdeltay) + 1)
            #
            if (iix>iib) and ((iix+1)<iie) and (iiy>ijb) and ((iiy+1)<ije):
               zxk = (zxk-dom2readgeom.geometry.ij2xy(iix,iiy)[0]*100000.) / zdeltax
               zyk = (zyk-dom2readgeom.geometry.ij2xy(iix,iiy)[1]*100000.) / zdeltay
               
               zvarcyl[jphi,jr]=zvarcart[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                zvarcart[iiy,iix+1]*zxk*(1-zyk)   +  \
                                zvarcart[iiy+1,iix]*(1-zxk)*zyk   +  \
                                zvarcart[iiy+1,iix+1]*zxk*zyk


    return zvarcyl

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def polar_calcwind(dom2readgeom,zucart,zvcart,iie,ije,icen,jcen,iradmax0,nphil,zdeltax,zdeltay,zdeltar,zdphi,hemisphere):

    ## *POLAR_CALC* - wind interpolation onto cylindrical grid
    ##
    ## PURPOSE
    ## -------
    # This subroutine interpolates horizontal wind U, V onto a cylindrical grid centred on (icen,jcen) at a given level (R=1 correspond to the TC's center). Here, one can get tangential Vt and radial winds, respectively.

    # INITIALIZATION
    #
    import numpy
    import epygram

    zradius=numpy.zeros(iradmax0, dtype=int)
    zucyl = numpy.zeros((nphil,iradmax0), dtype=float)
    zvcyl = numpy.zeros((nphil,iradmax0), dtype=float)
    zvtcyl = numpy.zeros((nphil,iradmax0), dtype=float)
    zurcyl = numpy.zeros((nphil,iradmax0), dtype=float)
    zuurcyl = numpy.zeros((nphil,iradmax0), dtype=float)
    zvurcyl = numpy.zeros((nphil,iradmax0), dtype=float)

    iib=0
    ijb=0


    for jr in range(1,iradmax0):
        zradius[jr] = jr * zdeltar
    #
    zxi0=dom2readgeom.geometry.ij2xy(icen,jcen)[0]*100000.
    zyj0=dom2readgeom.geometry.ij2xy(icen,jcen)[1]*100000.
    zx00=dom2readgeom.geometry.ij2xy(1,jcen)[0]*100000.
    zy00=dom2readgeom.geometry.ij2xy(icen,1)[1]*100000.

    for jr in range(iradmax0):
        for jphi in range(nphil):
            zphi = (jphi) * zdphi
            zxk  = zradius[jr] * numpy.cos(zphi) + zxi0
            zyk  = zradius[jr] * numpy.sin(zphi) + zyj0
            iix  = int(((zxk-zx00) / zdeltax) + 1 )
            iiy  = int(((zyk-zy00) / zdeltay) + 1)
            #
            if (iix>iib) and ((iix+1)<iie) and (iiy>ijb) and ((iiy+1)<ije):
               zxk = (zxk-dom2readgeom.geometry.ij2xy(iix,iiy)[0]*100000.) / zdeltax
               zyk = (zyk-dom2readgeom.geometry.ij2xy(iix,iiy)[1]*100000.) / zdeltay

               zucyl[jphi,jr]=zucart.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                              zucart.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                              zucart.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                              zucart.data[iiy+1,iix+1]*zxk*zyk
               zvcyl[jphi,jr]=zvcart.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                              zvcart.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                              zvcart.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                              zvcart.data[iiy+1,iix+1]*zxk*zyk              

               if hemisphere=='S':
                  zucyl[jphi,jr] =  -zucyl[jphi,jr]
                  zvcyl[jphi,jr] =  -zvcyl[jphi,jr]
               
               zvtcyl[jphi,jr] = zvcyl[jphi,jr]*numpy.cos(zphi) \
                                               - zucyl[jphi,jr]*numpy.sin(zphi)
               zurcyl[jphi,jr] = zucyl[jphi,jr]*numpy.cos(zphi) \
                                               + zvcyl[jphi,jr]*numpy.sin(zphi)
               
               if hemisphere=='S':
                  zurcyl[jphi,jr] = -zurcyl[jphi,jr]
               
               if zvtcyl[jphi,jr]<0:
                  zvtcyl[jphi,jr] = 0.               
            else:
                zvtcyl[jphi,jr]=-999.
                zurcyl[jphi,jr]=-999.

    return zvtcyl, zurcyl
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def mean_azim(zvarcyl,iradmax0,nphil,flag=-999.):

    import numpy 

    zvarcylMEAN = numpy.zeros(iradmax0, dtype=float)

    for jr in range(iradmax0):
        icount=0
        for jphi in range(0,nphil):
            if(zvarcyl[jphi,jr] != flag):
              zvarcylMEAN[jr] = zvarcylMEAN[jr] + zvarcyl[jphi,jr]
              icount = icount + 1
            else:
              zvarcylMEAN[jr]=flag
        if icount == len(range(0,nphil)):
           zvarcylMEAN[jr] = zvarcylMEAN[jr] / icount
        else:
           zvarcylMEAN[jr] = flag

    return zvarcylMEAN
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_RMW(zutcylmean,iradmax0,zdeltar):

    import numpy

    zrmwmax=0.

    for jr in range(iradmax0):
        if zutcylmean[jr]>zrmwmax:
           zrmwmax=zutcylmean[jr]
           irmw=jr
           zrmwMEAN=irmw*zdeltar / 1000.
    
    return zrmwMEAN
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_Vtmax_azim(zutcyl,iradmax0,zdeltar,nphil,zdphi):

    import numpy

    maxazim=numpy.max(zutcyl[:,:])

    return maxazim    
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_windshear(zullcylmean,zvllcylmean,zuulcylmean,zvulcylmean,iradmax0,zdeltar):

    import numpy

    zradius=numpy.zeros(iradmax0, dtype=int)

    for jr in range(iradmax0):
        zradius[jr] = jr * zdeltar


    ################################################################
    ### Compute the cylindrical area-weighted average of U and V ###
    ### Hanley et al, 2001, Corbosiero and Molteni 2002          ###
    ################################################################
    sumA=0.
    ncount=0
    uenvll=0.
    venvll=0.
    uenvul=0.
    venvul=0.
    if(zullcylmean[iradmax0-1]!=-999.):
      for jr in range(1,iradmax0):
          ncount=ncount=1
          A=numpy.math.pi*(zradius[jr]**2.-zradius[jr-1]**2.)
          sumA=sumA + A
          B=(zullcylmean[jr-1]+zullcylmean[jr])/2.
          C=(zvllcylmean[jr-1]+zvllcylmean[jr])/2.
          D=(zuulcylmean[jr-1]+zuulcylmean[jr])/2.
          E=(zvulcylmean[jr-1]+zvulcylmean[jr])/2.
          uenvll=uenvll+ B*A
          venvll=venvll+ C*A
          uenvul=uenvul+ D*A
          venvul=venvul+ E*A
    else:
      print("Attention à xradguess on sort du domaine")
      ffshear=-999.
      ddshear=-999.
    #  
    if(ncount!=0):
      uenvll=uenvll/sumA
      venvll=venvll/sumA
      uenvul=uenvul/sumA
      venvul=venvul/sumA


    Ushear=uenvul-uenvll
    Vshear=venvul-venvll  

    ffshear = numpy.sqrt(Ushear**2. + Vshear**2.)
    ffshear = ffshear * 1.944   # Conversion en knots
    ddshear = numpy.math.atan2(-Ushear,-Vshear)
    if(ddshear<0.):
      ddshear=ddshear+2*numpy.math.pi

    ddshear = ddshear * (180./numpy.math.pi)  

    return ffshear, ddshear 
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_THETAE(P,T,qv,iie,ije):

    import numpy

    zthetae=numpy.zeros((ije,iie), dtype=float)

    Lv=2400. # Coef de chaleur latente d'évaporation
    Cp=1004. # chaleur spécifique de l'air sec à pression constante
    R=287.   # constante spécifique de l'air sec dans le cadre de la loi des gaz parfaits
    K= R/Cp  # 

    for j in range(ije):
        for i in range(iie):
            zthetae[j,i] = (T.data[j,i] + (Lv/Cp)*qv[j,i]) * (1000./P.data[j,i])**K

    return zthetae
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_qv(P,T,Hu,iie,ije):

    import numpy

    zesl=numpy.zeros((ije,iie), dtype=float)
    zesi=numpy.zeros((ije,iie), dtype=float)
    ze=numpy.zeros((ije,iie), dtype=float)
    zqv=numpy.zeros((ije,iie), dtype=float)



    rv=461.5249 ; ra=287.0596                                          # Magnus-Tetens-Murray(67)-water & ice

    for j in range(ije):
        for i in range(iie):
            zesl[j,i]=610.78*numpy.math.exp( 17.269*(T.data[j,i]-273.15)/(T.data[j,i]-35.86) )          #water sat coeffs. units: P[Pa] T[K] zes[Pa]
            zesi[j,i]=610.78*numpy.math.exp( 21.875*(T.data[j,i]-273.15)/(T.data[j,i]-7.66 ) )                #ice sat coeffs. units: P[Pa] T[K] zes[Pa]

            ze[j,i] = Hu.data[j,i]/100.*zesl[j,i]                                                  # Hu(0-100) -> vapour partial pressure(Pa)

            zqv[j,i] = ze[j,i]/P.data[j,i]*ra/rv                                                   # zqv(kg/kg)  P(Pa)

    return zqv
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
