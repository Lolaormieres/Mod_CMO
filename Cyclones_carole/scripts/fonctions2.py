#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# import
from math import *
import numpy as np
import pandas as pd
import datetime
from Magics.macro import *

from scipy import arccos



######################################
# Transformations de formats de dates
#########################################

def transfo_date_collerToib(date, heure) :
   date, heure=str(date), str(heure)
   date2=date[0:4]+"-"+date[4:6]+"-"+date[6:8]+" "+heure+":00:00"
   return date2

def transfo_date_collerToseparer(date) :
   annee=str(date[0:4])
   mois=str(date[4:6])
   jour=str(date[6:8])
   return annee, mois, jour

def transfo_date_ibTosepare(date) :
   annee=str(date[0:4])
   mois=str(date[5:7])
   jour=str(date[8:10])
   heure=str(date[11:13])
   return annee, mois, jour, heure

def transfo_date_separeToib(annee, mois, jour, heure) :
   date=annee+"-"+mois+"-"+jour+" "+heure+":00:00"
   return date


######################################
# Calculs de dates
#########################################

def calcul_datefinale(datedebut, heure):

   datedeb=transfo_date_collerToib(datedebut, heure)
   datedeb0=np.datetime64(datedeb)
   datefin0=datedeb0+np.timedelta64(120, 'h')
   datefin1=np.datetime_as_string(datefin0)
   datefin=datefin1[:10]+' '+datefin1[11:]
   return datedeb, datefin



def calcul_date_a_une_Echeance(annee,mois,jour,heure, ech):

   datedeb=transfo_date_separeToib(annee,mois,jour, heure)
   datedeb0=np.datetime64(datedeb)
   datefin0=datedeb0+np.timedelta64(ech, 'h')
   datefin=np.datetime_as_string(datefin0)
   annee_fin, mois_fin, jour_fin, heure_fin = transfo_date_ibTosepare(datefin)
   return annee_fin, mois_fin, jour_fin, heure_fin


def calcul_date_a_une_Echeance_ib(annee,mois,jour,heure, ech):

   datedeb=transfo_date_separeToib(annee,mois,jour, heure)
   datedeb0=np.datetime64(datedeb)
   datefin0=datedeb0+np.timedelta64(ech, 'h')
   datefin1=np.datetime_as_string(datefin0)
   datefin=datefin1[:10]+' '+datefin1[11:]
   return datefin



#############################
# Récupération des données
#############################


# Dans ibtracks
def recup_data_ibtrack_between(datedeb, datefin):
   data0=np.loadtxt("ibtracs.last3years.list.v04r00.csv", dtype="str", delimiter=',', skiprows=2)
   datedeb=transfo_date_collerToib(datedeb)
   datefin=transfo_date_collerToib(datefin)
   data=data0[(datedeb <= data0[:,6]) & (data0[:,6] < datefin),:]
   return data



def recup_latlon_ibtrack_parTC(annee,mois, jour, heure, tcname) :
   datedeb=transfo_date_separeToib(annee, mois, jour, heure)
   data0=np.loadtxt("ibtracs.last3years.list.v04r00.csv", dtype="str", delimiter=',', skiprows=2)
   data_guess=data0[(data0[:,5]==tcname) & (data0[:,6]==datedeb), :]
   latguess=np.asarray(data_guess[:,8],dtype='float')
   longuess=np.asarray(data_guess[:,9],dtype='float')
   return latguess, longuess, 



#def recup_size_ibtrack_parTC(annee,mois, jour, heure) :
#   datedeb=transfo_date_separeToib(annee, mois, jour, heure)
#   data0=np.loadtxt("ibtracs.last3years.list.v04r00.csv", dtype="str", delimiter=',', skiprows=2)
#   data_guess=data0[(data0[:,27]==USA_34_NE) & (data0[:,28]==USA_34_NW) & (data0[:,29]==USA_34_SW) & (data0[:,30]==USA_34_SE) & (data0[:,6]==datedeb), :]
#   USA_34_NEguess=np.asarray(data_guess[:,14],dtype='float')
#   USA_34_NWguess=np.asarray(data_guess[:,15],dtype='float')
#   USA_34_SWguess=np.asarray(data_guess[:,16],dtype='float')
#   USA_34_SEguess=np.asarray(data_guess[:,17],dtype='float')
#   return tcname, USA_34_NEguess, USA_34_NWguess, USA_34_SWguess, USA_34_SEguess #, USA_64_NE, USA_64_NW, USA_64_SW, USA_64_SE


def recup_data_ibtrack(annee,mois, jour, heure, ibbasin, tcname) :
   datedeb=transfo_date_separeToib(annee, mois, jour, heure)
   data0=np.loadtxt("ibtracs.last3years.list.v04r00.csv", dtype="str", delimiter=',', skiprows=2)
   data_guess=data0[(data0[:,3]==ibbasin) & (data0[:,5]==tcname) & (data0[:,6]==datedeb), :]
   return data_guess





# Dans les fichiers de paramètres NAME_date.csv
def read_param_TC(tcname,annee_run, mois_run, jour_run, heure_run):
   data=np.loadtxt(tcname+'_'+annee_run+mois_run+jour_run+'_'+heure_run+'.csv',dtype="str", delimiter=',')
   return data


def latlon_param_TC_parmembre_pour1ech(tcname,annee_run, mois_run, jour_run, heure_run, ech, membre) :
   data0=read_param_TC(tcname,annee_run, mois_run, jour_run, heure_run)
   data=data0[(data0[:,1] == membre) & (data0[:,2] == str(ech)), :]
   lat=np.asarray(data[:,3], dtype='float')
   lon=np.asarray(data[:,4], dtype='float')
   return lat, lon


def vmaxpmin_param_TC_pour1ech(tcname,annee_run, mois_run, jour_run, heure_run, ech) :
   data0=read_param_TC(tcname,annee_run, mois_run, jour_run, heure_run)
   data=data0[(data0[:,2] == str(ech)), :]
   vmax=np.asarray(data[:,5], dtype='float')
   pmin=np.asarray(data[:,6], dtype='float')
   return vmax, pmin
#plus bon maintenant

#def exstension_param_TC_pour1ech(tcname,annee_run, mois_run, jour_run, heure_run, ech) :
#   data0=read_param_TC(tcname,annee_run, mois_run, jour_run, heure_run)
#   data=data0[(data0[:,2] == str(ech)), :]
#   zrmw34obsNE=np.asarray(data[:,7], dtype='float')
#   zrmw34obsNW=np.asarray(data[:,8], dtype='float')
#   zrmw34obsSW=np.asarray(data[:,9], dtype='float')
#   zrmw34obsSE=np.asarray(data[:,10], dtype='float')
#   return zrmw34obsNE, zrmw34obsNW, zrmw34obsSW, zrmw34obsSE

def equiv_param_fichiername(nom_du_fichier):

  tcname=nom_du_fichier[:-16]
  annee_run=nom_du_fichier[-15:-11]
  mois_run=nom_du_fichier[-11:-9]
  jour_run=nom_du_fichier[-9:-7]
  heure_run=nom_du_fichier[-6:-4]

  return tcname, annee_run, mois_run, jour_run, heure_run


def equiv_colparam(param) :
  if param == "lat" : 
    num_col = 3
    num_col_ib=8
  elif param == "lon" : 
    num_col = 4
    num_col_ib=9
  elif param == "vmax850" : 
    num_col = 5
    num_col_ib=23
  elif param == "vmax925" :
    num_col = 6
    num_col_ib=23
  elif param == "vmax10m" : 
    num_col = 7
    num_col_ib=23
  elif param == "vmax50m" : 
    num_col = 8 
    num_col_ib=23
  elif param == "vmax100m" : 
    num_col = 9 
    num_col_ib=23 
  elif param == "vmax500m" : 
    num_col = 10 
    num_col_ib=23
  elif param == "vmaxgust1h" : 
    num_col = 11 
    num_col_ib=23
  elif param == "vmaxgust3h" : 
    num_col = 12 
    num_col_ib=23
  elif param == "pmin" : 
    num_col = 13
    num_col_ib=24
  else : 
    print("pas de paramètre : "+param)
  return num_col, num_col_ib




# Récupération des noms de fichiers

def recup_fichiers_dans1repertoire(rep) :

   from os import listdir
   fichiers=[]
   for f in listdir(rep) :
      if ".csv" in f : fichiers.append(f)
   print(fichiers)
   return fichiers

 
####################
# calcul nbr_tc et leur noms pour 1 date

def nombre_TC(annee, mois, jour, heure, ibbasin):
   date=transfo_date_separeToib(annee, mois, jour, heure)
   data=np.loadtxt("ibtracs.last3years.list.v04r00.csv", dtype="str", delimiter=',', skiprows=2)
   nbr_tc=0
   list_tc=[]
   for jnb in range(len(data)) :
      if data[jnb,6] == date and data[jnb,3] == ibbasin : 
         nbr_tc += 1
         list_tc.append(data[jnb,5])

   print('Il y a '+str(nbr_tc)+' cyclones analysés sur le bassin '+ibbasin+' pour le '+jour+'/'+mois+'/'+annee+' à '+heure+' UTC')   
   return nbr_tc, list_tc





#####################################################
# Fonction qui fait correspondre analyse et prévision
######################################################

def correspondance_analyse_prevision(tcname,annee_run, mois_run, jour_run, heure_run, ech, membre) :
   annee, mois, jour, heure = calcul_date_a_une_Echeance(annee_run,mois_run,jour_run, heure_run, ech)
   lat_prevision, lon_prevision = latlon_param_TC_parmembre_pour1ech(tcname,annee_run, mois_run, jour_run, heure_run, ech, membre)
   lat_analyse, lon_analyse= recup_latlon_ibtrack_parTC(annee,mois, jour, heure, tcname)
   return lat_analyse, lon_analyse, lat_prevision , lon_prevision
   





#############################
# plots, affectation couleur, taille, ...

def classification_TC_parmembre(data, member, ech):
   vmax, pmin, =[],[]
   colorTC, HeightTC, ThickTC, styleTC=[],[],[],[]
   vmax[0:ech]=data.loc[data.loc[:,1] == member,5]
   pmin[0:ech]=data.loc[data.loc[:,1] == member,6]

   for jc in range(ech) :
      symbTC=15
      if (pmin[jc]>=1010. or vmax[jc]<30.):                            # Non classifié
         colorTC.append("black")
         HeightTC.append(0.1)
         ThickTC.append(4)
         styleTC.append("dash")
      else :
         HeightTC.append(0.1)
         ThickTC.append(4)
         styleTC.append("solid")
         if (pmin[jc]<1010 and vmax[jc]>=136.):                           # H5          
            colorTC.append("RGB(128,0,128)" )
         if (pmin[jc]<1010 and vmax[jc]>=114. and vmax[jc]<136.):         # H4
            colorTC.append("RGB(220,20,60)")
         if (pmin[jc]<1010 and vmax[jc]>=96. and vmax[jc]<114):           # H3
            colorTC.append("RGB(255,69,0)")
         if (pmin[jc]<1010 and vmax[jc]>=83. and vmax[jc]<96):            # H2
            colorTC.append("RGB(255,140,0)")
         if (pmin[jc]<1010 and vmax[jc]>=64. and vmax[jc]<83):            # H1
            colorTC.append("RGB(255,255,0)")  
         if (pmin[jc]<1010 and vmax[jc]>=35. and vmax[jc]<64):            # TS
            colorTC.append("RGB(173,255,47)")
         if (pmin[jc]<1010 and vmax[jc]>=30. and vmax[jc]<35):            # TD
            colorTC.append("RGB(30,144,255)")   

   return colorTC, HeightTC, ThickTC, styleTC, symbTC





# choix des membres

def choix_des_membres_pour1ech(tcname,annee_run, mois_run, jour_run, heure_run, ech):

   liste_complete=['mb000', 'mb001', 'mb002', 'mb003', 'mb004', 'mb005', 'mb006', 'mb007', 'mb008', 'mb009', 'mb010','mb011', 'mb012', 'mb013', 'mb014','mb015', 'mb016', 'mb017','mb018', 'mb019', 'mb020', 'mb021', 'mb022', 'mb023', 'mb024', 'mb025', 'mb026', 'mb027', 'mb028', 'mb029', 'mb030', 'mb031', 'mb032', 'mb033', 'mb034']
   liste_membres_valides=[]
   vmax, pmin = vmaxpmin_param_TC_pour1ech(tcname,annee_run, mois_run, jour_run, heure_run, ech)

   for jc in range(len(vmax)) :
      if (pmin[jc]<1010. and vmax[jc]>=30.) :
         liste_membres_valides.append(liste_complete[jc])

   if len(liste_membres_valides) < 8 : liste_membres_valides=[]

   return liste_membres_valides
#ok






# taille boite, position legende equivalent bassin dans ibtracks en fonction du bassin choisi

def attributs_bassin(ibbassin) :
   attributs =dict()

   if ibbassin=='NA' :
     attributs["Hemisphere"] = 'N' 
     attributs["basin"] = "ATLN"
     attributs["latn"]=44.
     attributs["lats"]=4.
     attributs["lono"]=-90.
     attributs["lone"]=-14.
     attributs["pos_legend"] = minput(
                       input_x_values    =    [-17.5],
                       input_y_values    =    [40.5]
                             )

   elif ibbassin=='EP':
     attributs["Hemisphere"] = 'N'
     attributs["basin"] = "PCNE"
     attributs["latn"]=44.
     attributs["lats"]=4.
     attributs["lono"]=-180.
     attributs["lone"]=-80.
     attributs["pos_legend"] = minput(
                       input_x_values    =    [-85.],
                       input_y_values    =    [40.5]
                             ) 

   elif ibbassin=='WP' :
     attributs["Hemisphere"] = 'N'
     attributs["basin"] = "PCNO"
     attributs["latn"]=44.
     attributs["lats"]=4.
     attributs["lono"]=100.
     attributs["lone"]=180. 
     attributs["pos_legend"] = minput(
                       input_x_values    =    [104.],
                       input_y_values    =    [40.5]
                            )

   elif ibbassin=='NI':
     attributs["Hemisphere"] = 'N' 
     attributs["basin"] = "NOI"
     attributs["latn"]=44.
     attributs["lats"]=4.
     attributs["lono"]=40.
     attributs["lone"]=100.
     attributs["pos_legend"] = minput(
                       input_x_values    =    [42.5],
                       input_y_values    =    [40.5]
                            )  
   
   elif ibbassin=='SI':
     attributs["Hemisphere"] = 'S'
     attributs["basin"] = "SEOI ou SOOI"
     attributs["latn"]=0.
     attributs["lats"]=-40.
     attributs["lono"]=30.
     attributs["lone"]=95.
     attributs["pos_legend"] = minput(
                       input_x_values    =    [88.5],
                       input_y_values    =    [-36.5]
                             )

   elif ibbassin=='SP':
     attributs["Hemisphere"] = 'S'
     attributs["basin"] = "PCS"
     attributs["latn"]=0.
     attributs["lats"]=-40.
     attributs["lono"]=146.92
     attributs["lone"]=180.
     attributs["pos_legend"] = minput(
                       input_x_values    =    [150.],
                       input_y_values    =    [-35.]
                             )  

   return attributs
# ok


'''
elif ibbassin=='SI':
     attributs["Hemisphere"] = 'S'
     attributs["basin"] = "SEOI"
     attributs["latn"]=0.
     attributs["lats"]=-40.
     attributs["lono"]=95.
     attributs["lone"]=146.92 
     attributs["pos_legend"] = minput(
                       input_x_values    =    [143.5],
                       input_y_values    =    [-25.5]
                             )
'''














