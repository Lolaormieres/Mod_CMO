#!/usr/bin/env python
#!/bin/env python

from Magics.macro import *
import epygram
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

epygram.init_env()



#####################################################################
# Initializations
#####################################################################

fd=open('Hurr-Track-Aro.nam','r')
data_nam=numpy.loadtxt('Hurr-Track-Aro.nam',dtype='str')[:]

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
echmax=int(data_nam[12])
expe=data_nam[13]

LGRIB2=True

zdeltax = 2500.
zdeltay = 2500.


print("On se trouve sur le bassin=", bassin)

if(domaine=='GLOB025'):
  zdeltax = 25000.
  zdeltay = 25000.
  iie=1440
  ije=721

if(domaine=='CARAIB0025'):
  iie=945
  ije=529

if(domaine=='INDIEN0025'):
  iie=1395
  ije=747

if(domaine=='atllarge0025'):
  iie=1201
  ije=841

if(domaine=='mdr0025'):
  iie=1761
  ije=441

if(domaine=='nclarge0025'):
  iie=1201
  ije=521

list_ech=[]
echstep=6
for t in range(0,echmax+echstep,echstep):
    list_ech.append(t) # Les niveaux lus dans les gribs Arome
ech=len(list_ech)

print('On a '+str(ech)+' echeances a traiter')
list_tc=[]


#list_ech=[3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48]
#list_ech=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48]
#ech=len(list_ech)
#print('On a '+str(ech)+' echeances a traiter')
#list_tc=[]


Alltrack=[]
Alldate=[]
Allshear=[]
AnaPos=[]
AllTrace=[]
All64wind=[]
All34wind=[]
All50wind=[]
#Allmsg108=[]
#Allmsg108no=[]
Allwind10=[]
AllTimepos=[]
AllObs=[]

speedguess = 50.
xradguess = 200. # en km
nphil=361  # number of azimuts (for cyclonic projection)
nthel=9   # number of azimuts (for drawing of wind quadrants)
xboxwind=30. # 
zspeedmin=100.
iib=0
ijb=0

latinc=5
loninc=5


zdeltar = max(zdeltax,zdeltay)
zdphi = numpy.math.pi / 16.
zdtheta = numpy.math.pi / 16.

iradmax0 = numpy.int64((xradguess*1000. / zdeltar) + 2)



### For Visualisation
radmax     = xradguess
windmax = 200. # (en kts)
rmwmax = xradguess #(en kms)
minpres = 900
maxpres = 1020

#####################################################################
# Paramètres pour le domaine de tracé
#####################################################################
area = mmap(subpage_map_projection="cylindrical",
        subpage_lower_left_longitude=lono,
        subpage_lower_left_latitude=lats,
        subpage_upper_right_longitude=lone,
        subpage_upper_right_latitude=latn,
         ) 
full_screen = mmap(
 #                   subpage_x_position       = 2.,
 #                   subpage_y_position       = 2.,
                 subpage_x_length         = 26.,
 #                   subpage_y_length         = 10.,
                 page_x_position          = 0.,
                 page_y_position          = 0.,
                 page_x_length            = 40.,
                 page_y_length            = 20.,
                 page_id_line             = 'off'
           )

coast = mcoast(map_coastline_land_shade = "on",
           #map_coastline_land_shade_colour = "RGB(191,191,191)",
           map_coastline_land_shade_colour = "RGB(255,248,229)",
           map_coastline_sea_shade = "on",
           #map_coastline_sea_shade_colour = "RGB(100,162,202)",
           map_coastline_sea_shade_colour = "RGB(166,191,221)",
           map_coastline_thickness = 3,
           map_grid_line_style = "dot",
           map_grid_thickness = 1,
           map_grid_colour = "black",
           map_label = "on",
           map_label_height = 0.9,
           map_grid_latitude_increment = latinc,
           map_grid_longitude_increment = loninc,
           map_coastline_colour = "black")


coast_noshade = mcoast(map_coastline_land_shade = "on",
           #map_coastline_land_shade_colour = "RGB(191,191,191)",
           map_coastline_land_shade_colour = "RGB(255,248,229)",
           map_coastline_sea_shade = "off",
           #map_coastline_sea_shade_colour = "RGB(100,162,202)",
           map_coastline_sea_shade_colour = "RGB(166,191,221)",
           map_coastline_thickness = 3,
           map_grid_line_style = "dot",
           map_grid_thickness = 1,
           map_grid_colour = "black",
           map_label = "on",
           map_label_height = 0.9,
           map_grid_latitude_increment = latinc,
           map_grid_longitude_increment = loninc,
           map_coastline_colour = "black")




if(bassin=='CUSTOM-HN'):
   Hemisphere = 'N' 
   pos_legend1 = minput(
                       input_x_values    =    [-96.5],
                       input_y_values    =    [16.]
                             )
   symb_legend1 = msymb(
                         symbol_image_format = "png",
	                     symbol_image_path = "legendewindmax_bw.png",
	                     symbol_type = "marker",
	                     symbol_marker_mode = "image",
	                     symbol_image_width = 7,
                         symbol_image_height = 2,
	                     )
   loc_1 = minput(
                       input_x_values    =    [-61.],
                       input_y_values    =    [16.2]
                             )
   symb_loc1 = msymb(
                       symbol_type    =    "both",
                       symbol_text_list    =    ["Guadeloupe"],
                       symbol_marker_index    =    2,
                       symbol_colour    =    "black",
                       symbol_height    =    0.,
                       symbol_text_font_size = 0.8,
                       symbol_text_font_colour = "black",
                       symbol_text_position = "right",
                       symbol_text_font_style = "italic",
                       )
   loc_2 = minput(
                       input_x_values    =    [-60.5],
                       input_y_values    =    [14.6]
                             )
   symb_loc2 = msymb(
                       symbol_type    =    "both",
                       symbol_text_list    =    ["Martinique"],
                       symbol_marker_index    =    2,
                       symbol_colour    =    "black",
                       symbol_height    =    0.,
                       symbol_text_font_size = 0.8,
                       symbol_text_font_colour = "black",
                       symbol_text_position = "right",
                       symbol_text_font_style = "italic",
                       )
   loc_3 = minput(
                       input_x_values    =    [-62.7],
                       input_y_values    =    [18.1]
                             )
   symb_loc3 = msymb(
                       symbol_type    =    "both",
                       symbol_text_list    =    ["St Barth & St Martin"],
                       symbol_marker_index    =    2,
                       symbol_colour    =    "black",
                       symbol_height    =    0.,
                       symbol_text_font_size = 0.8,
                       symbol_text_font_colour = "black",
                       symbol_text_position = "right",
                       symbol_text_font_style = "italic",
                       )

if(bassin=='CUSTOM-HS'):
   Hemisphere = 'S'
   pos_legend1 = minput(
                       input_x_values    =    [172.],
                       input_y_values    =    [-19.5]
                             )
   symb_legend1 = msymb(
                         symbol_image_format = "png",
                             symbol_image_path = "legendewindmax_bw.png",
                             symbol_type = "marker",
                             symbol_marker_mode = "image",
                             symbol_image_width = 4,
                         symbol_image_height = 1,
                             )
   loc_1 = minput(
                       input_x_values    =    [56.],
                       input_y_values    =    [-21.5]
                             )
   symb_loc1 = msymb(
                       symbol_type    =    "both",
                       symbol_text_list    =    ["La Réunion"],
                       symbol_marker_index    =    2,
                       symbol_colour    =    "black",
                       symbol_height    =    0.,
                       symbol_text_font_size = 0.8,
                       symbol_text_font_colour = "black",
                       symbol_text_position = "right",
                       symbol_text_font_style = "italic",
                       )
   loc_2 = minput(
                       input_x_values    =    [58.],
                       input_y_values    =    [-20.]
                             )
   symb_loc2 = msymb(
                       symbol_type    =    "both",
                       symbol_text_list    =    ["Maurice"],
                       symbol_marker_index    =    2,
                       symbol_colour    =    "black",
                       symbol_height    =    0.,
                       symbol_text_font_size = 0.8,
                       symbol_text_font_colour = "black",
                       symbol_text_position = "right",
                       symbol_text_font_style = "italic",
                       )
   loc_3 = minput(
                       input_x_values    =    [45.3],
                       input_y_values    =    [-13.]
                             )
   symb_loc3 = msymb(
                       symbol_type    =    "both",
                       symbol_text_list    =    ["Mayotte"],
                       symbol_marker_index    =    2,
                       symbol_colour    =    "black",
                       symbol_height    =    0.,
                       symbol_text_font_size = 0.8,
                       symbol_text_font_colour = "black",
                       symbol_text_position = "right",
                       symbol_text_font_style = "italic",
                       )


#AllTrace.append(loc_1)
#AllTrace.append(symb_loc1)
#AllTrace.append(loc_2)
#AllTrace.append(symb_loc2)
#AllTrace.append(loc_3)
#AllTrace.append(symb_loc3)
AllTrace.append(pos_legend1)
AllTrace.append(symb_legend1)


title = mtext(text_lines=['Trajectoire - Intensité ARPEGE: '+annee+'-'+mois+'-'+jour+' '+heure+':00 UTC'],
                        text_font_size = 1.5,
                        text_colour= "charcoal",
                        text_font_style='bold',
                )



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

    zucyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zvcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zmod10cyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zprepcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zprep1cyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zprep2cyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)
    zutcyl = numpy.zeros((ech,nphil,iradmax0), dtype=float)

    zutcylMEANAZIM = numpy.zeros((ech,iradmax0), dtype=float)
    zutcylNE = numpy.zeros((ech,iradmax0), dtype=float)
    zutcylNW = numpy.zeros((ech,iradmax0), dtype=float)
    zutcylSW = numpy.zeros((ech,iradmax0), dtype=float)
    zutcylSE = numpy.zeros((ech,iradmax0), dtype=float)
    zprecylMEANAZIM = numpy.zeros((ech,iradmax0), dtype=float)
    zu850cylMEANAZIM = numpy.zeros((ech,iradmax0), dtype=float)
    zv850cylMEANAZIM = numpy.zeros((ech,iradmax0), dtype=float)

    zmod10cylMAX = numpy.zeros((ech,iradmax0), dtype=float)
    zwindMAX = numpy.zeros((ech), dtype=float)
    zpmin = numpy.zeros((ech), dtype=float)
    zrmwMEAN = numpy.zeros((ech), dtype=float)
    zrmwMEANTS = numpy.zeros((ech), dtype=float)
    zrmwMEANHU = numpy.zeros((ech), dtype=float)
    zrmwMEANMHU = numpy.zeros((ech), dtype=float)
    #
    ## Write the data file for scores
    if os.path.exists(model+'-'+annee+mois+jour+heure+'-'+tcname):
       os.remove(model+'-'+annee+mois+jour+heure+'-'+tcname)
    fout=open(model+'-'+annee+mois+jour+heure+'-'+tcname,'a')
#
    # Begin of the loop on time jc  
    for jc in range(ech):
        zspeedmin=100.
        zrmwmax=0.
        zrmwmaxNE=0.
        zrmwmaxNW=0.
        zrmwmaxSW=0.
        zrmwmaxSE=0.
        #lecture  
        f = epygram.formats.resource('grid.'+model+'-forecast.'+domaine+'+'+str(list_ech[jc]).zfill(4)+':00.grib', 'r', 'GRIB')
        if(jc-1>=0):
          fprec = epygram.formats.resource('grid.'+model+'-forecast.'+domaine+'+'+str(list_ech[jc-1]).zfill(4)+':00.grib', 'r', 'GRIB')  
        if (LGRIB2==True):
           mslp = f.readfield({'name': 'Pressure reduced to MSL'})
           u10=f.readfield({'name':'10 metre U wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
           v10=f.readfield({'name':'10 metre V wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
           mod10=numpy.sqrt(u10.data**2. + v10.data**2.)
           if(jc-1>=0):
             rain2 = f.readfield({'parameterCategory': 1, 'parameterNumber':76})
             convrain = f.readfield({'parameterCategory': 1, 'parameterNumber':76})
             lsrain = f.readfield({'parameterCategory': 1, 'parameterNumber':77})
             rain2.data[:,:] = convrain.data[:,:] + lsrain.data[:,:]
             if(list_ech[jc]==echstep):
                rain1 = fprec.readfield({'name': 'Pressure reduced to MSL'})
                rain1.data[:,:] = 0.
             else:
                rain1 = fprec.readfield({'parameterCategory': 1, 'parameterNumber':76})
                convrain = fprec.readfield({'parameterCategory': 1, 'parameterNumber':76})
                lsrain = fprec.readfield({'parameterCategory': 1, 'parameterNumber':77})
                rain1.data[:,:] = convrain.data[:,:] + lsrain.data[:,:]
                #rain1 = fprec.readfield({'name': 'Rain precipitation rate'})
           else:
             if(list_ech[jc]>0.):
               rain2 = f.readfield({'parameterCategory': 1, 'parameterNumber':76})
               convrain = f.readfield({'parameterCategory': 1, 'parameterNumber':76})
               lsrain = f.readfield({'parameterCategory': 1, 'parameterNumber':77})
               rain2.data[:,:] = convrain.data[:,:] + lsrain.data[:,:]
               #rain2 = f.readfield({'name': 'Rain precipitation rate'})
             else:
               rain2 = f.readfield({'name': 'Pressure reduced to MSL'})
               rain2.data[:,:] = 0.
        else:
           mslp = f.readfield({'indicatorOfParameter': 2})
           if(jc-1>=0):
             rain2 = f.readfield({'indicatorOfParameter': 150}) 
             rain1 = fprec.readfield({'indicatorOfParameter': 150})
           else:
             rain2 = f.readfield({'indicatorOfParameter': 150})
        mslp.operation('/',100.)
        ztime_basis[0]=mslp.validity.getbasis()
        ztime_term=mslp.validity.term()
        ztime[jc]=mslp.validity[0].get()
        #
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
        ##########################################################################
        # Détermination de la zone de caclul des conditions cycloniques 
        # (on considère un cyclone potentiel de déplaçant à une vitesse speedguess
        ##########################################################################
        if(jc-1>=0):
           data=numpy.loadtxt('pos.geo'+str(list_ech[jc-1]).zfill(2)+'-'+tcname)[:]
           latpos_1=data[0]
           lonpos_1=data[1]
        else:
           latguess=latobs
           longuess=lonobs 
           latpos_1=latguess
           lonpos_1=longuess
           if(bassin=='CUSTOM-HN'):
             print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % latguess)+' °N et '+str("% 4.2f" % (longuess*-1.))+' °W')
           if(bassin=='CUSTOM-HS'):
             print('Le cyclone '+tcname+' est analysé par '+str("% 4.2f" % (latguess*-1.))+' °S et '+str("% 4.2f" % longuess)+' °E')
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
           centre_ll=mslp.geometry.ij2ll(iicen,ijcen)
#   
           #####################################################################
           # Détermination du centre cyclonique iicen et ijcen 
           # (minimum de vent dans l'oeil)
           #####################################################################
#
           if (LGRIB2==True):
              u10=f.readfield({'name':'10 metre U wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
              v10=f.readfield({'name':'10 metre V wind component', 'typeOfLevel':'heightAboveGround', 'level': 10})
           else:
              u10=f.readfield({'indicatorOfTypeOfLevel': 105, 'level': 10, 'table2Version': 1, 'editionNumber': 1, 'indicatorOfParameter': 33})
              v10=f.readfield({'indicatorOfTypeOfLevel': 105, 'level': 10, 'table2Version': 1, 'editionNumber': 1, 'indicatorOfParameter': 34})

           mod10=numpy.sqrt(u10.data**2. + v10.data**2.)   
              
#
           iavgwindi = int(round(xboxwind*1000./zdeltax))
           iavgwindj = int(round(xboxwind*1000./zdeltay))
#
           for jj in range(max(ijcen-iavgwindj,1),min(ijcen+iavgwindj,ije)):
               for ji in range(max(iicen-iavgwindi,1),min(iicen+iavgwindi,iie)):
                   zspeed=numpy.sqrt(u10.data[jj,ji]**2. + v10.data[jj,ji]**2.)
                   if zspeed<= zspeedmin:
                      zspeedmin=zspeed
                      iicen2 = ji
                      ijcen2 = jj
#  
           iicen=iicen2
           ijcen=ijcen2
#
           centre_ll=mslp.geometry.ij2ll(iicen,ijcen)
           if(bassin=='CUSTOM-HN'):
              print('A l\'échéance '+str(list_ech[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % centre_ll[1])+' °N et'+str("% 5.2f" % ((centre_ll[0]-360.)*-1)+' °W'))
           if(bassin=='CUSTOM-HS'):
              print('A l\'échéance '+str(list_ech[jc]).zfill(2)+' '+tcname+' est relocalisé par'+str("% 7.2f" % (centre_ll[1]*-1.))+' °S et'+str("% 5.2f" % (centre_ll[0]-360.)+' °E'))   
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
                      zucyl[jc,jphi,jr]=u10.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                        u10.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                        u10.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                        u10.data[iiy+1,iix+1]*zxk*zyk
                      zvcyl[jc,jphi,jr]=v10.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                        v10.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                        v10.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                        v10.data[iiy+1,iix+1]*zxk*zyk
                      zmod10cyl[jc,jphi,jr]=mod10.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                        mod10.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                        mod10.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                        mod10.data[iiy+1,iix+1]*zxk*zyk                  
                      if(jc-1>=0):                  
                         zprep1cyl[jc,jphi,jr]=rain1.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                               rain1.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                               rain1.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                               rain1.data[iiy+1,iix+1]*zxk*zyk 
                         zprep2cyl[jc,jphi,jr]=rain2.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                               rain2.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                               rain2.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                               rain2.data[iiy+1,iix+1]*zxk*zyk
                         zprepcyl[jc,jphi,jr] = zprep2cyl[jc,jphi,jr] - zprep1cyl[jc,jphi,jr]  
                      else:
                         zprep2cyl[jc,jphi,jr]=rain2.data[iiy,iix]*(1-zxk)*(1-zyk) +  \
                                               rain2.data[iiy,iix+1]*zxk*(1-zyk)   +  \
                                               rain2.data[iiy+1,iix]*(1-zxk)*zyk   +  \
                                               rain2.data[iiy+1,iix+1]*zxk*zyk
                         zprepcyl[jc,jphi,jr] = zprep2cyl[jc,jphi,jr]                        

                      if Hemisphere=='S':
                         zucyl[jc,jphi,jr] =  -zucyl[jc,jphi,jr]
                         zvcyl[jc,jphi,jr] =  -zvcyl[jc,jphi,jr]
                      zutcyl[jc,jphi,jr] = zvcyl[jc,jphi,jr]*numpy.cos(zphi) \
                                           - zucyl[jc,jphi,jr]*numpy.sin(zphi)
                                                   
                      if zutcyl[jc,jphi,jr]<0: 
                         zutcyl[jc,jphi,jr] = 0.
                      zutcyl[jc,jphi,jr] = zutcyl[jc,jphi,jr] * 1.944      # Conversion en kts
                      zucyl[jc,jphi,jr] = zucyl[jc,jphi,jr] * 1.944      # Conversion en kts
                      zvcyl[jc,jphi,jr] = zvcyl[jc,jphi,jr] * 1.944      # Conversion en kts
                      zmod10cyl[jc,jphi,jr] = zmod10cyl[jc,jphi,jr] * 1.944      # Conversion en kts
                   
                   #
                   else:
                      zutcyl[jc,jphi,jr]=-999. 
                      zucyl[jc,jphi,jr]=-999.
                      zvcyl[jc,jphi,jr]=-999.
                      zprepcyl[jc,jphi,jr]=-999.
               icount=0
               zmod10cylMAX[jc,jr] = numpy.max(zmod10cyl[jc,:,jr])
            #
            #
               icount=0
               for jphi in range(0,361):
                   if(zutcyl[jc,jphi,jr]!=-999.):
                      zutcylMEANAZIM[jc,jr] = zutcylMEANAZIM[jc,jr] + zutcyl[jc,jphi,jr]
                      zprecylMEANAZIM[jc,jr] = zprecylMEANAZIM[jc,jr] + zprepcyl[jc,jphi,jr]
                      icount = icount + 1
                   else:
                      zutcylMEANAZIM[jc,jr]=-999.
                      zprecylMEANAZIM[jc,jr]=-999.
               if icount == len(range(0,361)):
                  zutcylMEANAZIM[jc,jr] = zutcylMEANAZIM[jc,jr] / icount
                  zprecylMEANAZIM[jc,jr] = zprecylMEANAZIM[jc,jr] / icount
               else:
                  zutcylMEANAZIM[jc,jr] = -999.
                  zprecylMEANAZIM[jc,jr] = -999.
            #
            #
               icount=0
               for jphi in range(0,91):
                   if(zutcyl[jc,jphi,jr]!=-999.):
                      zutcylNE[jc,jr] = zutcylNE[jc,jr] + zutcyl[jc,jphi,jr]
                      icount = icount + 1
                   else:
                      zutcylNE[jc,jr]=-999.
               if icount == len(range(0,91)):
                  zutcylNE[jc,jr] = zutcylNE[jc,jr] / icount
               else:
                  zutcylNE[jc,jr] = -999.
            #
            #
               icount=0
               for jphi in range(90,181):
                   if(zutcyl[jc,jphi,jr]!=-999.):
                      zutcylNW[jc,jr] = zutcylNW[jc,jr] + zutcyl[jc,jphi,jr]
                      icount = icount + 1
                   else:
                      zutcylNW[jc,jr]=-999.
               if icount == len(range(90,181)):
                  zutcylNW[jc,jr] = zutcylNW[jc,jr] / icount
               else:
                  zutcylNW[jc,jr] = -999.            
            #
            #
               icount=0
               for jphi in range(180,271):
                   if(zutcyl[jc,jphi,jr]!=-999.):
                      zutcylSW[jc,jr] = zutcylSW[jc,jr] + zutcyl[jc,jphi,jr]
                      icount = icount + 1
                   else:
                      zutcylSW[jc,jr]=-999.
               if icount == len(range(180,271)):
                  zutcylSW[jc,jr] = zutcylSW[jc,jr] / icount
               else:
                  zutcylSW[jc,jr] = -999.
            #
            #
               icount=0
               for jphi in range(270,361):
                   if(zutcyl[jc,jphi,jr]!=-999.):
                      zutcylSE[jc,jr] = zutcylSE[jc,jr] + zutcyl[jc,jphi,jr]
                      icount = icount + 1
                   else:
                      zutcylSE[jc,jr]=-999.
               if icount == len(range(270,361)):
                  zutcylSE[jc,jr] = zutcylSE[jc,jr] / icount
               else:
                  zutcylSE[jc,jr] = -999.            
            
#        
           Ushear = 0.
           Vshear = 0.
           Modshear = 0.
#        
           if os.path.exists('track.csv-'+str(list_ech[jc]).zfill(2)+'-'+tcname):
              os.remove('track.csv-'+str(list_ech[jc]).zfill(2)+'-'+tcname)
           f=open('track.csv-'+str(list_ech[jc]).zfill(2)+'-'+tcname,'w')
           lat=centre_ll[1]
           lon=centre_ll[0]-360.
           #
           #
           date='date(2012-10-23)'
           vmax=numpy.max(zmod10cylMAX[jc,:])
           zwindMAX[jc]=vmax 
           pmin=mslpmin
           zpmin[jc]=pmin
           #
           irmw=0
           irmwNE=0
           irmwNW=0
           irmwSW=0
           irmwSE=0
           lat34NE=[lat]
           lon34NE=[lon]
           lat34NW=[lat]
           lon34NW=[lon]
           lat34SW=[lat]
           lon34SW=[lon]
           lat34SE=[lat]
           lon34SE=[lon]
           lat50NE=[lat]
           lon50NE=[lon]
           lat50NW=[lat]
           lon50NW=[lon]
           lat50SW=[lat]
           lon50SW=[lon]
           lat50SE=[lat]
           lon50SE=[lon]
           lat64NE=[lat]
           lon64NE=[lon]
           lat64NW=[lat]
           lon64NW=[lon]
           lat64SW=[lat]
           lon64SW=[lon]
           lat64SE=[lat]
           lon64SE=[lon]
           for jr in range(iradmax0):
              if zutcylMEANAZIM[jc,jr]>zrmwmax:
                zrmwmax=zutcylMEANAZIM[jc,jr]
                irmw=jr
                zrmwMEAN[jc]=irmw*zdeltar / 1000.
           #
           for jr in range(iradmax0):
              if zutcylNE[jc,jr]>zrmwmaxNE:
                zrmwmaxNE=zutcylNE[jc,jr]
                irmwNE=jr
           #
           for jr in range(iradmax0):
              if zutcylNW[jc,jr]>zrmwmaxNW:
                zrmwmaxNW=zutcylNW[jc,jr]
                irmwNW=jr
           for jr in range(iradmax0):
              if zutcylSW[jc,jr]>zrmwmaxSW:
                zrmwmaxSW=zutcylSW[jc,jr]
                irmwSW=jr
           #
           for jr in range(iradmax0):
              if zutcylSE[jc,jr]>zrmwmaxSE:
                zrmwmaxSE=zutcylSE[jc,jr]
                irmwSE=jr
           #     
           zrmwMEANTS[jc]=-999.
           zrmw50NE=-999.
           zrmw64NE=-999.
           zrmw34NE=-999.
           zrmw50NW=-999.
           zrmw64NW=-999.
           zrmw34NW=-999.
           zrmw50SW=-999.
           zrmw64SW=-999.
           zrmw34SW=-999.
           zrmw50SE=-999.
           zrmw64SE=-999.
           zrmw34SE=-999.
           if (irmw>0):
              if (zutcylMEANAZIM[jc,irmw]>34):
                for jr in range(irmw,iradmax0):
                   if (zutcylMEANAZIM[jc,jr]<34 and zutcylMEANAZIM[jc,jr]!=-999.):
                     irmwts=jr
                     zrmwMEANTS[jc]=irmw*zdeltar / 1000. + (irmwts-irmw)*zdeltar / 1000.
                     break
           if (irmwNE>0):
              if (zutcylNE[jc,irmwNE]>50):
                for jr in range(irmwNE,iradmax0):
                   if (zutcylNE[jc,jr]<50 and zutcylNE[jc,jr]!=-999.):
                     irmw50=jr
                     zrmw50NE=irmwNE*zdeltar / 1000. + (irmw50-irmwNE)*zdeltar / 1000.
                     lat50NE=[]
                     lon50NE=[]
                     for jtheta in range(0,9):
                        ztheta = (jtheta) * zdtheta
                        lat50NE.append(lat+(zrmw50NE*numpy.sin(ztheta)/111.))
                        lon50NE.append(lon+(zrmw50NE*numpy.cos(ztheta)/111.))
                     break
           if (irmwNE>0):
              if (zutcylNE[jc,irmwNE]>64):
                for jr in range(irmwNE,iradmax0):
                   if (zutcylNE[jc,jr]<64 and zutcylNE[jc,jr]!=-999.):
                     irmw64=jr
                     zrmw64NE=irmwNE*zdeltar / 1000. + (irmw64-irmwNE)*zdeltar / 1000.
                     lat64NE=[]
                     lon64NE=[]
                     for jtheta in range(0,9):
                        ztheta = (jtheta) * zdtheta
                        lat64NE.append(lat+(zrmw64NE*numpy.sin(ztheta)/111.))
                        lon64NE.append(lon+(zrmw64NE*numpy.cos(ztheta)/111.))
                     break
           if (irmwNE>0):
              if (zutcylNE[jc,irmwNE]>34):
                for jr in range(irmwNE,iradmax0):
                   if (zutcylNE[jc,jr]<34 and zutcylNE[jc,jr]!=-999.):
                     irmw34=jr
                     zrmw34NE=irmwNE*zdeltar / 1000. + (irmw34-irmwNE)*zdeltar / 1000.
                     lat34NE=[]
                     lon34NE=[]
                     for jtheta in range(0,9):
                        ztheta = (jtheta) * zdtheta
                        lat34NE.append(lat+(zrmw34NE*numpy.sin(ztheta)/111.))
                        lon34NE.append(lon+(zrmw34NE*numpy.cos(ztheta)/111.))
                     break                 
           #
           if (irmwNW>0):
              if (zutcylNW[jc,irmwNW]>50):
                for jr in range(irmwNW,iradmax0):
                   if (zutcylNW[jc,jr]<50 and zutcylNW[jc,jr]!=-999.):
                     irmw50=jr
                     zrmw50NW=irmwNW*zdeltar / 1000. + (irmw50-irmwNW)*zdeltar / 1000.
                     lat50NW=[]
                     lon50NW=[]
                     for jtheta in range(8,17):
                        ztheta = (jtheta) * zdtheta
                        lat50NW.append(lat+(zrmw50NW*numpy.sin(ztheta)/111.))
                        lon50NW.append(lon+(zrmw50NW*numpy.cos(ztheta)/111.))
                     break
           if (irmwNW>0):
              if (zutcylNW[jc,irmwNW]>64):
                for jr in range(irmwNW,iradmax0):
                   if (zutcylNW[jc,jr]<64 and zutcylNW[jc,jr]!=-999.):
                     irmw64=jr
                     zrmw64NW=irmwNW*zdeltar / 1000. + (irmw64-irmwNW)*zdeltar / 1000.
                     lat64NW=[]
                     lon64NW=[]
                     for jtheta in range(8,17):
                        ztheta = (jtheta) * zdtheta
                        lat64NW.append(lat+(zrmw64NW*numpy.sin(ztheta)/111.))
                        lon64NW.append(lon+(zrmw64NW*numpy.cos(ztheta)/111.))
                     break
           if (irmwNW>0):
              if (zutcylNW[jc,irmwNW]>34):
                for jr in range(irmwNW,iradmax0):
                   if (zutcylNW[jc,jr]<34 and zutcylNW[jc,jr]!=-999.):
                     irmw34=jr
                     zrmw34NW=irmwNW*zdeltar / 1000. + (irmw34-irmwNW)*zdeltar / 1000.
                     lat34NW=[]
                     lon34NW=[]
                     for jtheta in range(8,17):
                        ztheta = (jtheta) * zdtheta
                        lat34NW.append(lat+(zrmw34NW*numpy.sin(ztheta)/111.))
                        lon34NW.append(lon+(zrmw34NW*numpy.cos(ztheta)/111.))
                     break
           #
           if (irmwSW>0):
              if (zutcylSW[jc,irmwSW]>50):
                for jr in range(irmwSW,iradmax0):
                   if (zutcylSW[jc,jr]<50 and zutcylSW[jc,jr]!=-999.):
                     irmw50=jr
                     zrmw50SW=irmwSW*zdeltar / 1000. + (irmw50-irmwSW)*zdeltar / 1000.
                     lat50SW=[]
                     lon50SW=[]
                     for jtheta in range(16,25):
                        ztheta = (jtheta) * zdtheta
                        lat50SW.append(lat+(zrmw50SW*numpy.sin(ztheta)/111.))
                        lon50SW.append(lon+(zrmw50SW*numpy.cos(ztheta)/111.))
                     break
           if (irmwSW>0):
              if (zutcylSW[jc,irmwSW]>64):
                for jr in range(irmwSW,iradmax0):
                   if (zutcylSW[jc,jr]<64 and zutcylSW[jc,jr]!=-999.):
                     irmw64=jr
                     zrmw64SW=irmwSW*zdeltar / 1000. + (irmw64-irmwSW)*zdeltar / 1000.
                     lat64SW=[]
                     lon64SW=[]
                     for jtheta in range(16,25):
                        ztheta = (jtheta) * zdtheta
                        lat64SW.append(lat+(zrmw64SW*numpy.sin(ztheta)/111.))
                        lon64SW.append(lon+(zrmw64SW*numpy.cos(ztheta)/111.))
                     break
           if (irmwSW>0):
              if (zutcylSW[jc,irmwSW]>34):
                for jr in range(irmwSW,iradmax0):
                   if (zutcylSW[jc,jr]<34 and zutcylSW[jc,jr]!=-999.):
                     irmw34=jr
                     zrmw34SW=irmwSW*zdeltar / 1000. + (irmw34-irmwSW)*zdeltar / 1000.
                     lat34SW=[]
                     lon34SW=[]
                     for jtheta in range(16,25):
                        ztheta = (jtheta) * zdtheta
                        lat34SW.append(lat+(zrmw34SW*numpy.sin(ztheta)/111.))
                        lon34SW.append(lon+(zrmw34SW*numpy.cos(ztheta)/111.))
                     break 
           #
           if (irmwSE>0):
              if (zutcylSE[jc,irmwSE]>50):
                for jr in range(irmwSE,iradmax0):
                   if (zutcylSE[jc,jr]<50 and zutcylSE[jc,jr]!=-999.):
                     irmw50=jr
                     zrmw50SE=irmwSE*zdeltar / 1000. + (irmw50-irmwSE)*zdeltar / 1000.
                     lat50SE=[]
                     lon50SE=[]
                     for jtheta in range(24,33):
                        ztheta = (jtheta) * zdtheta
                        lat50SE.append(lat+(zrmw50SE*numpy.sin(ztheta)/111.))
                        lon50SE.append(lon+(zrmw50SE*numpy.cos(ztheta)/111.))
                     break
           if (irmwSE>0):
              if (zutcylSE[jc,irmwSE]>64):
                for jr in range(irmwSE,iradmax0):
                   if (zutcylSE[jc,jr]<64 and zutcylSE[jc,jr]!=-999.):
                     irmw64=jr
                     zrmw64SE=irmwSE*zdeltar / 1000. + (irmw64-irmwSE)*zdeltar / 1000.
                     lat64SE=[]
                     lon64SE=[]
                     for jtheta in range(24,33):
                        ztheta = (jtheta) * zdtheta
                        lat64SE.append(lat+(zrmw64SE*numpy.sin(ztheta)/111.))
                        lon64SE.append(lon+(zrmw64SE*numpy.cos(ztheta)/111.))
                     break
           if (irmwSE>0):
              if (zutcylSE[jc,irmwSE]>34):
                for jr in range(irmwSE,iradmax0):
                   if (zutcylSE[jc,jr]<34 and zutcylSE[jc,jr]!=-999.):
                     irmw34=jr
                     zrmw34SE=irmwSE*zdeltar / 1000. + (irmw34-irmwSE)*zdeltar / 1000.
                     lat34SE=[]
                     lon34SE=[]
                     for jtheta in range(24,33):
                        ztheta = (jtheta) * zdtheta
                        lat34SE.append(lat+(zrmw34SE*numpy.sin(ztheta)/111.))
                        lon34SE.append(lon+(zrmw34SE*numpy.cos(ztheta)/111.))
                     break    
           #
           if(zrmw34NE!=-999. and zrmw34NW!=-999.):
             lat34NE.append(lat+(zrmw34NW*numpy.sin(numpy.math.pi/2.)/111.))
             lon34NE.append(lon+(zrmw34NW*numpy.cos(numpy.math.pi/2.)/111.))
           if(zrmw34NW!=-999. and zrmw34SW!=-999.):
             lat34NW.append(lat+(zrmw34SW*numpy.sin(numpy.math.pi)/111.))
             lon34NW.append(lon+(zrmw34SW*numpy.cos(numpy.math.pi)/111.)) 
           if(zrmw34SW!=-999. and zrmw34SE!=-999.):
             lat34SW.append(lat+(zrmw34SE*numpy.sin(3.*numpy.math.pi/2.)/111.))
             lon34SW.append(lon+(zrmw34SE*numpy.cos(3.*numpy.math.pi/2.)/111.))  
           if(zrmw34SE!=-999. and zrmw34NE!=-999.):
             lat34SE.append(lat+(zrmw34NE*numpy.sin(0.)/111.))
             lon34SE.append(lon+(zrmw34NE*numpy.cos(0.)/111.)) 
           #
           #
           if(zrmw50NE!=-999. and zrmw50NW!=-999.):
             lat50NE.append(lat+(zrmw50NW*numpy.sin(numpy.math.pi/2.)/111.))
             lon50NE.append(lon+(zrmw50NW*numpy.cos(numpy.math.pi/2.)/111.))
           if(zrmw50NW!=-999. and zrmw50SW!=-999.):
             lat50NW.append(lat+(zrmw50SW*numpy.sin(numpy.math.pi)/111.))
             lon50NW.append(lon+(zrmw50SW*numpy.cos(numpy.math.pi)/111.)) 
           if(zrmw50SW!=-999. and zrmw50SE!=-999.):
             lat50SW.append(lat+(zrmw50SE*numpy.sin(3.*numpy.math.pi/2.)/111.))
             lon50SW.append(lon+(zrmw50SE*numpy.cos(3.*numpy.math.pi/2.)/111.))  
           if(zrmw50SE!=-999. and zrmw50NE!=-999.):
             lat50SE.append(lat+(zrmw50NE*numpy.sin(0.)/111.))
             lon50SE.append(lon+(zrmw50NE*numpy.cos(0.)/111.))
           #
           #
           if(zrmw64NE!=-999. and zrmw64NW!=-999.):
             lat64NE.append(lat+(zrmw64NW*numpy.sin(numpy.math.pi/2.)/111.))
             lon64NE.append(lon+(zrmw64NW*numpy.cos(numpy.math.pi/2.)/111.))
           if(zrmw64NW!=-999. and zrmw64SW!=-999.):
             lat64NW.append(lat+(zrmw64SW*numpy.sin(numpy.math.pi)/111.))
             lon64NW.append(lon+(zrmw64SW*numpy.cos(numpy.math.pi)/111.))
           if(zrmw64SW!=-999. and zrmw64SE!=-999.):
             lat64SW.append(lat+(zrmw64SE*numpy.sin(3.*numpy.math.pi/2.)/111.))
             lon64SW.append(lon+(zrmw64SE*numpy.cos(3.*numpy.math.pi/2.)/111.))
           if(zrmw64SE!=-999. and zrmw64NE!=-999.):
             lat64SE.append(lat+(zrmw64NE*numpy.sin(0.)/111.))
             lon64SE.append(lon+(zrmw64NE*numpy.cos(0.)/111.))
           #
           if(zrmw34SE==-999.):
             lat34SW.append(lat)
             lon34SW.append(lon)
             lat34NE.insert(0,lat)
             lon34NE.insert(0,lon)
           if(zrmw34SW==-999.):
             lat34NW.append(lat)
             lon34NW.append(lon)
             lat34SE.insert(0,lat)
             lon34SE.insert(0,lon)
           if(zrmw34NW==-999.):
             lat34NE.append(lat)
             lon34NE.append(lon)
             lat34SW.insert(0,lat)
             lon34SW.insert(0,lon)
           if(zrmw34NE==-999.):
             lat34SE.append(lat)
             lon34SE.append(lon)
             lat34NW.insert(0,lat)
             lon34NW.insert(0,lon)
           #
           if(zrmw50SE==-999.):
             lat50SW.append(lat)
             lon50SW.append(lon)
             lat50NE.insert(0,lat)
             lon50NE.insert(0,lon)
           if(zrmw50SW==-999.):
             lat50NW.append(lat)
             lon50NW.append(lon)
             lat50SE.insert(0,lat)
             lon50SE.insert(0,lon)
           if(zrmw50NW==-999.):
             lat50NE.append(lat)
             lon50NE.append(lon)
             lat50SW.insert(0,lat)
             lon50SW.insert(0,lon)
           if(zrmw50NE==-999.):
             lat50SE.append(lat)
             lon50SE.append(lon)
             lat50NW.insert(0,lat)
             lon50NW.insert(0,lon)
           #  
           if(zrmw64SE==-999.):
             lat64SW.append(lat)
             lon64SW.append(lon)
             lat64NE.insert(0,lat)
             lon64NE.insert(0,lon)
           if(zrmw64SW==-999.):
             lat64NW.append(lat)
             lon64NW.append(lon)
             lat64SE.insert(0,lat)
             lon64SE.insert(0,lon)
           if(zrmw64NW==-999.):
             lat64NE.append(lat)
             lon64NE.append(lon)
             lat64SW.insert(0,lat)
             lon64SW.insert(0,lon)
           if(zrmw64NE==-999.):
             lat64SE.append(lat)
             lon64SE.append(lon)
             lat64NW.insert(0,lat)
             lon64NW.insert(0,lon)

           ####################
           #
           zrmwMEANHU[jc]=-999.
           if (irmw>0):
              if (zutcylMEANAZIM[jc,irmw]>64):
                for jr in range(irmw,iradmax0):
                   if (zutcylMEANAZIM[jc,jr]<64 and zutcylMEANAZIM[jc,jr]!=-999.):
                     irmwhu=jr 
                     zrmwMEANHU[jc]=irmw*zdeltar / 1000. + (irmwhu-irmw)*zdeltar / 1000.
                     break
           zrmwMEANMHU[jc]=-999.
           if (irmw>0):
              if (zutcylMEANAZIM[jc,irmw]>96):
                for jr in range(irmw,iradmax0):
                   if (zutcylMEANAZIM[jc,jr]<96 and zutcylMEANAZIM[jc,jr]!=-999.):
                     irmwmhu=jr
                     zrmwMEANMHU[jc]=irmw*zdeltar / 1000. + (irmwmhu-irmw)*zdeltar / 1000.
                     break  
           #print('A l\'échéance '+str(list_ech[jc]).zfill(2)+' Extension HU= '+str(zrmwMEANHU[jc]).zfill(2))
           #f.write(str("%4.2f" % lat)+','+str("%4.2f" % lon)+','+date+','+str("%4i" % heure)+','+str("%2i" % vmax)+','+str("%4i" % pmin )+',')
           f.write(str("%4.2f" % lat)+','+str("%4.2f" % lon)+','+date+','+heure+','+str("%2i" % vmax)+','+str("%4i" % pmin )+',')
           f.close()
#
           if os.path.exists('pos.geo'+str(list_ech[jc]).zfill(2)+'-'+tcname):
              os.remove('pos.geo'+str(list_ech[jc]).zfill(2)+'-'+tcname)
           f=open('pos.geo'+str(list_ech[jc]).zfill(2)+'-'+tcname,'a')
           f.write("#GEO\n")
           f.write("#FORMAT XY_VECTOR\n")
           f.write("# lat lon height date time u v\n")
           f.write("#DATA\n")
           f.write(str("%4.2f" % lat)+' '+str("%4.2f" % lon)+' 0 20030617 1200 '+str("%4.2f" % Ushear)+' '+str("%4.2f" % Vshear)+"\n")
           f.close()
#   
           if(bassin=='CUSTOM-HN'):
             if(lat>latn or lat<lats or lon<lono or lon>lone):
                zutcylMEANAZIM[jc,:] = -999.
                zprecylMEANAZIM[jc,:] = -999.


           if(bassin=='CUSTOM-HS'):
             if(lat>latn or lat<lats or lon<lono or lon>lone):
                zutcylMEANAZIM[jc,:] = -999.
                zprecylMEANAZIM[jc,:] = -999.
           #print("Valeur de iradmax0=",iradmax0)
           #print("Valeur de irmw=",irmw)
           jrmin = int(((rmwmax / 3.) / zdeltar*1000.) + 1) 
           #if(zutcylMEANAZIM[jc,jrmin]==-999.):
           #  zutcylMEANAZIM[jc,:]=-999.
           #  zprecylMEANAZIM[jc,:] = -999.
           #  break  # 
#
           ThickTC=18
           ThickTCW=4
           ThickTCW96=4
           styleTC="solid"
           HeightTC=0.8
           symbTC=18
           if (vmax>=137.):
              symbTC=15                              # H5
              colorTC="RGBA(128,0,128,1.)"
           if (vmax>=113. and vmax<137.):
              symbTC=15                               # H4
              colorTC="RGBA(220,20,60,1.)"
           if (vmax>=96. and vmax<113):
              symbTC=15                               # H3
              colorTC="RGBA(255,69,0,1.)"
           if (vmax>=83. and vmax<96):
              symbTC=15                                # H2
              colorTC="RGBA(255,140,0,1.)"
           if (vmax>=64. and vmax<83):
              symbTC=15                                # H1
              colorTC="RGBA(255,255,0,1.)"
           if (vmax>=34. and vmax<64):
              symbTC=15                                # TS
              colorTC="RGBA(173,255,47,1.)"
           if (vmax>=30. and vmax<34):
              symbTC=15                                # TD
              colorTC="RGBA(30,144,255,1.)"
           if (vmax<30.):
              symbTC=15
              colorTC="RGBA(0,0,0,0.4)"
           
           if Modshear<=15.*0.514:
              colorFlag="yellowish_green"
           if (Modshear>15.*0.514 and Modshear<25.*0.514):
              colorFlag="orange_yellow"   
           if Modshear>=25.*0.514:
              colorFlag="red_orange"
#         
           if(jc-1>=0):
              data = mtable(table_filename = "track.csv-"+str(list_ech[jc-1]).zfill(2)+'-'+tcname,
                      table_variable_identifier_type='index',
                      table_latitude_variable = "1",
                      table_longitude_variable = "2",
                      table_value_variable = "-1",
                      table_header_row = 0,
                      )
           if(list_ech[jc]==72):
              data = mtable(table_filename = "track.csv-"+str(list_ech[jc]).zfill(2)+'-'+tcname,
                      table_variable_identifier_type='index',
                      table_latitude_variable = "1",
                      table_longitude_variable = "2",
                      table_value_variable = "-1",
                      table_header_row = 0,
                      )
           if(list_ech[jc]==0):
              data = mtable(table_filename = "track.csv-"+str(list_ech[jc]).zfill(2)+'-'+tcname,
                      table_variable_identifier_type='index',
                      table_latitude_variable = "1",
                      table_longitude_variable = "2",
                      table_value_variable = "-1",
                      table_header_row = 0,
                      )   
           tc_pos = minput(
                          input_x_values    =    [lon],
                          input_y_values    =    [lat]
                                )
           symb_tc_pos = msymb(
                         symbol_type    =    "both",
                         symbol_marker_index    =    15,
                         symbol_colour    =    colorTC,
                         symbol_height    =    0.,
                         symbol_text_font_size = 0.,
                         )
           vent = mgeo(geo_input_file_name = "pos.geo"+str(list_ech[jc]).zfill(2)+'-'+tcname)
          
           if(jc-1>=0):
              track = minput(
                       input_x_values    =    [lonpos_1,lon],
                       input_y_values    =    [latpos_1,lat]
                       )
           else:
              track = minput(
                       input_x_values    =    [lon,lon],
                       input_y_values    =    [lat,lat]
                       )

           trackline = mgraph(
                graph_line_colour    =    colorTC,
                graph_line_thickness    =    ThickTC,
                graph_line_style    =   styleTC
                )
           track64NE = minput(
                       input_x_values    =    lon64NE,
                       input_y_values    =    lat64NE
                       )
           track50NE = minput(
                       input_x_values    =    lon50NE,
                       input_y_values    =    lat50NE
                       )
           track34NE = minput(
                       input_x_values    =    lon34NE,
                       input_y_values    =    lat34NE
                       )
           trackline64NE = mgraph(
                graph_line_colour    =    "RGBA(255,255,0,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           trackline50NE = mgraph(
                graph_line_colour    =    "RGBA(143,31,63,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           trackline34NE = mgraph(
                graph_line_colour    =    "RGBA(173,255,47,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           track64NW = minput(
                       input_x_values    =    lon64NW,
                       input_y_values    =    lat64NW
                       )
           track50NW = minput(
                       input_x_values    =    lon50NW,
                       input_y_values    =    lat50NW
                       )
           track34NW = minput(
                       input_x_values    =    lon34NW,
                       input_y_values    =    lat34NW
                       )
           trackline64NW = mgraph(
                graph_line_colour    =    "RGBA(255,255,0,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           trackline50NW = mgraph(
                graph_line_colour    =    "RGBA(143,31,63,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           trackline34NW = mgraph(
                graph_line_colour    =    "RGBA(173,255,47,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )  
           track64SW = minput(
                       input_x_values    =    lon64SW,
                       input_y_values    =    lat64SW
                       )
           track50SW = minput(
                       input_x_values    =    lon50SW,
                       input_y_values    =    lat50SW
                       )
           track34SW = minput(
                       input_x_values    =    lon34SW,
                       input_y_values    =    lat34SW
                       )
           trackline64SW = mgraph(
                graph_line_colour    =    "RGBA(255,255,0,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           trackline50SW = mgraph(
                graph_line_colour    =    "RGBA(143,31,63,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           trackline34SW = mgraph(
                graph_line_colour    =    "RGBA(173,255,47,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           track64SE = minput(
                       input_x_values    =    lon64SE,
                       input_y_values    =    lat64SE
                       )
           track50SE = minput(
                       input_x_values    =    lon50SE,
                       input_y_values    =    lat50SE
                       )
           track34SE = minput(
                       input_x_values    =    lon34SE,
                       input_y_values    =    lat34SE
                       )
           trackline64SE = mgraph(
                graph_line_colour    =    "RGBA(255,255,0,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           trackline50SE = mgraph(
                graph_line_colour    =    "RGBA(143,31,63,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )
           trackline34SE = mgraph(
                graph_line_colour    =    "RGBA(173,255,47,1.)",
                graph_line_thickness    =    ThickTCW,
                graph_line_style    =   styleTC
                )            

           line=msymb(
                   symbol_type="marker",
                   symbol_marker_index = symbTC,
                   symbol_colour = colorTC,
                   symbol_outline="on",
                   symbol_outline_thickness = 2,
                   symbol_height = 0.,
                   symbol_text_font_colour = "black",       
                   symbol_connect_line ="on",
                   symbol_text_font_size = 0.
                   )
           pos0h=msymb(
                       symbol_type="marker",
                       symbol_marker_index = 3,
                       symbol_colour = "RGB(0,0,0)",
                       symbol_outline="on",
                       symbol_outline_thickness = 2,
                       symbol_outline_colour = "black",
                       symbol_height = 0.,
                       symbol_text_font_colour = "black",
                       symbol_text_font_size = 0.8,
                       symbol_text_list    =  ["00"]
                       )
           pos12h=msymb(
                       symbol_type="marker",
                       symbol_marker_index = 3,
                       symbol_colour = "RGB(0,0,0)",
                       symbol_outline="on",
                       symbol_outline_thickness = 2,
                       symbol_outline_colour = "black",
                       symbol_height = 0.,
                       symbol_text_font_colour = "black",
                       symbol_text_font_size = 0.8,
                       symbol_text_list    =  ["12"]
                       )
           pos24h=msymb(
                       symbol_type="marker",
                       symbol_marker_index = 3,
                       symbol_colour = "RGB(0,0,0)",
                       symbol_outline="on",
                       symbol_outline_thickness = 2,
                       symbol_outline_colour = "black",
                       symbol_height = 0.,
                       symbol_text_font_colour = "red",
                       symbol_text_font_size = 1.,
                       symbol_text_list    =  ["24"]
                       )
           pos36h=msymb(
                       symbol_type="marker",
                       symbol_marker_index = 3,
                       symbol_colour = "RGB(0,0,0)",
                       symbol_outline="on",
                       symbol_outline_thickness = 2,
                       symbol_outline_colour = "black",
                       symbol_height = 0.,
                       symbol_text_font_colour = "black",
                       symbol_text_font_size = 0.8,
                       symbol_text_list    =  ["36"]
                       )
           pos48h=msymb(
                       symbol_type="marker",
                       symbol_marker_index = 3,
                       symbol_colour = "RGB(170,170,170)",
                       symbol_outline="on",
                       symbol_outline_thickness = 2,
                       symbol_outline_colour = "black",
                       symbol_height = 0.,
                       symbol_text_font_colour = "red",
                       symbol_text_font_size = 1.,
                       symbol_text_list    =  ["48"]
                       )
           pos60h=msymb(
                       symbol_type="marker",
                       symbol_marker_index = 3,
                       symbol_colour = "RGB(170,170,170)",
                       symbol_outline="on",
                       symbol_outline_thickness = 2,
                       symbol_outline_colour = "black",
                       symbol_height = 0.,
                       symbol_text_font_colour = "black",
                       symbol_text_font_size = 0.8,
                       symbol_text_list    =  ["60"]
                       )
           pos72h=msymb(
                       symbol_type="marker",
                       symbol_marker_index = 3,
                       symbol_colour = "RGB(170,170,170)",
                       symbol_outline="on",
                       symbol_outline_thickness = 2,
                       symbol_outline_colour = "black",
                       symbol_height = 0.,
                       symbol_text_font_colour = "red",
                       symbol_text_font_size = 1.,
                       symbol_text_list    =  ["72"]
                       )
           pos96h=msymb(
                       symbol_type="marker",
                       symbol_marker_index = 3,
                       symbol_colour = "RGB(170,170,170)",
                       symbol_outline="on",
                       symbol_outline_thickness = 2,
                       symbol_outline_colour = "black",
                       symbol_height = 0.,
                       symbol_text_font_colour = "red",
                       symbol_text_font_size = 1.,
                       symbol_text_list    =  ["96"]
                       )
           #Pointing the analyzed tc
           if(jc==0):
              tc_ana_pos = minput(
                          input_x_values    =    [longuess],
                          input_y_values    =    [latguess]
                                ) 
              #definition of the symbol
              if(bassin=='CUSTOM-HN'):
                 symb_ana = msymb(
                             symbol_type    =    "both",
                             symbol_text_list    =    [tcname],
                             symbol_marker_index    =    15,
                             symbol_colour    =    "navy",
                             symbol_height    =    0.,
                             symbol_text_font_size = 1.,
                             symbol_text_font_colour = "navy",
                             symbol_text_position = "bottom",
                             symbol_text_font_style = "bold",
                             )
              if(bassin=='CUSTOM-HS'):
                 symb_ana = msymb(
                             symbol_type    =    "both",
                             symbol_text_list    =    [tcname],
                             symbol_marker_index    =    18,
                             symbol_colour    =    "navy",
                             symbol_height    =    0.4,
                             symbol_text_font_size = 1.,
                             symbol_text_font_colour = "navy",
                             symbol_text_position = "top",
                             symbol_text_font_style = "italic",
                             )   
              AnaPos.append(tc_ana_pos)
              AnaPos.append(symb_ana)
           windshear=mwind(wind_field_type='flags',
                   wind_flag_colour= colorFlag,
                   wind_flag_thickness = 2,
                   wind_flag_length = 0.6,
                   wind_flag_origin_marker = "off"
                   )
           
           if(zutcylMEANAZIM[jc,jrmin]!=-999.):
              #if(numpy.max(zutcylMEANAZIM[jc,:])>20.):
              #print("Valeur de vent moyen azim max=",numpy.max(zutcylMEANAZIM[jc,:])) 
              Alltrack.append(track)
              Alltrack.append(trackline)
              if(((jc-1>=0 and (list_ech[jc]%12==0)) or (list_ech[jc]==0))):
                 Alltrack.append(data)
                 Alltrack.append(line)
              #Allshear.append(vent)
              #if(list_ech[jc]==0):
              #       AllTimepos.append(data)
              #       AllTimepos.append(pos0h)
              #if(list_ech[jc]==12):
              #       AllTimepos.append(data)
              #       AllTimepos.append(pos12h)
              if(list_ech[jc]==24):
                     AllTimepos.append(data)
                     AllTimepos.append(pos24h)
              #if(list_ech[jc]==36):
              #       AllTimepos.append(data)
              #       AllTimepos.append(pos36h)
              if(list_ech[jc]==48):
                     AllTimepos.append(data)
                     AllTimepos.append(pos48h)
              #if(list_ech[jc]==60):
              #       AllTimepos.append(data)
              #       AllTimepos.append(pos60h)       
              if(list_ech[jc]==72):
                     AllTimepos.append(data)
                     AllTimepos.append(pos72h)       
              if(list_ech[jc]==96):
                     AllTimepos.append(data)
                     AllTimepos.append(pos96h)       
              if(list_ech[jc]%12==0):
                 All34wind.append(track34NE)
                 All34wind.append(trackline34NE)
                 All34wind.append(track34NW)
                 All34wind.append(trackline34NW)
                 All34wind.append(track34SW)
                 All34wind.append(trackline34SW)
                 All34wind.append(track34SE)
                 All34wind.append(trackline34SE)
                 #All50wind.append(track50NE)
                 #All50wind.append(trackline50NE)
                 #All50wind.append(track50NW)
                 #All50wind.append(trackline50NW)
                 #All50wind.append(track50SW)
                 #All50wind.append(trackline50SW)
                 #All50wind.append(track50SE)
                 #All50wind.append(trackline50SE)
                 All64wind.append(track64NE)
                 All64wind.append(trackline64NE)
                 All64wind.append(track64NW)
                 All64wind.append(trackline64NW)
                 All64wind.append(track64SW)
                 All64wind.append(trackline64SW)
                 All64wind.append(track64SE)
                 All64wind.append(trackline64SE)
                #Allshear.append(windshear)
        else:
           zutcylMEANAZIM[jc,:]=-999.
           zprecylMEANAZIM[jc,:] = -999.
           zpmin[jc]=-999.
           zwindMAX[jc]=-999.
        #
        if(list_ech[jc]%2==0):
          fout.write(tcname+' '+str(list_ech[jc]).zfill(4)+' '+str(ztime[jc])+' '+str("%4.2f" % lat)+' '+str("%4.2f" % lon)+' '+str("%2i" % vmax)+' '+str("%4i" % pmin )+"\r\n")
        #fout.write(tcname+' '+str(ztime[jc])+' '+str("%4.2f" % lat)+' '+str("%4.2f" % lon)+"\r\n")
        #
        if(numpy.max(zutcylMEANAZIM[jc,:])==-999.):
           zpmin[jc]=-999.
           zwindMAX[jc]=-999.
           zrmwMEAN[jc]=-999.
           zrmwMEANTS[jc]=-999.
           zrmwMEANHU[jc]=-999.
           zrmwMEANMHU[jc]=-999.
        #print("Valeur de pmin",zpmin[jc])   
        # End of the loop on time jc
    
    #         
    ### Adding track and wind observations

    fobs=open('BT_HISTORY_'+tcname,'r')
    data_obs=numpy.loadtxt('BT_HISTORY_'+tcname,dtype='str')[:,:]


    nbr_obs = 0
    while fobs.readline():
        nbr_obs += 1
    print('Il y a '+str(nbr_obs)+' positions analysC)es pour le cyclone '+tcname)

    for jobs in range(nbr_obs-1):
        dateobs=data_obs[jobs,1]
        heureobs=data_obs[jobs,2]
        latobs=float(data_obs[jobs,3])
        lonobs=float(data_obs[jobs,4])
        windobs=float(data_obs[jobs,5]) 

        HeightTC=0.5
        if (windobs>=137.):
           symbTC=16                              # H5
           colorTC="RGB(128,0,128)"
        if (windobs>=113. and windobs<137.):
           symbTC=16                               # H4
           colorTC="RGB(220,20,60)"
        if (windobs>=96. and windobs<113):
           symbTC=16                               # H3
           colorTC="RGB(255,69,0)"
        if (windobs>=83. and windobs<96):
           symbTC=16                                # H2
           colorTC="RGB(255,140,0)"
        if (windobs>=64. and windobs<83):
           symbTC=16                                # H1
           colorTC="RGB(255,255,0)"
        if (windobs>=34. and windobs<64):
           symbTC=16                                # TS
           colorTC="RGB(173,255,47)"
        if (windobs>=30. and windobs<34):
           symbTC=16                                # TD
           colorTC="RGB(30,144,255)"
        if (windobs<30.):
           symbTC=16                              # Non classifiC)
           colorTC="black"
           ThickTC=4
           styleTC="dash"

        obs_pos = minput(
                           input_x_values    =    [lonobs],
                           input_y_values    =    [latobs]
                           )

        symb_obs=msymb(
                   symbol_type="marker",
                   symbol_marker_index = 18,
                   symbol_colour = colorTC,
                   symbol_outline="on",
                   symbol_outline_thickness = 3,
                   symbol_height = HeightTC,
                   symbol_text_font_colour = "black",
                   symbol_connect_line ="on",
                   symbol_text_font_size = 0.
                   )

        pos_date = minput(
                              input_x_values    =    [lonobs+0.5],
                              input_y_values    =    [latobs]
                                   )
        symb_date = msymb(
                           symbol_type    =    "both",
                           symbol_text_list    =    [dateobs+' '+heureobs],
                           symbol_marker_index    =    3,
                           symbol_colour    =    "black",
                           symbol_height    =    0.,
                           symbol_text_font_size = 0.6,
                           symbol_text_font_colour = "black",
                           symbol_text_position = "right",
                           symbol_text_font_style = "bold",
                           )
        if(jobs%4==0):
           Alldate.append(pos_date)
           Alldate.append(symb_date)
        #if(jobs%2==0):  
        AllObs.append(obs_pos)
        AllObs.append(symb_obs)   

fout.close()

#setting the output
output1 = output(   
                output_formats = ['png'],
                output_name = expe+"track-radiuswind-"+bassin+"-"+str(ztime_basis[0]),
                output_width = 1400,                                              
                output_name_first_page_number = "off"
        )

#output2 = output(
#                output_formats = ['png'],
#                output_name = "track-cloud-"+bassin+"-"+str(ztime_basis[0]),
#                output_width = 1400,
#                output_name_first_page_number = "off"
#        )




#plot(output1, area, full_screen, coast, All34wind, All50wind, All64wind, Alltrack, AnaPos, AllTrace, Alldate, Allshear, AllTimepos, title)
plot(output1, area, full_screen, coast, Alltrack, AllObs, AnaPos, AllTrace, Alldate, Allshear, AllTimepos, title)
#plot(output2, area, full_screen, coast, Allmsg108, Allwind10, Allmsg108no, coast_noshade, Alltrack, AnaPos, AllTrace, Alldate, Allshear, title)

