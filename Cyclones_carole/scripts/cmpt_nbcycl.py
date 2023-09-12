#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Import

import epygram
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os 

epygram.init_env()

from fonctions2 import nombre_TC, recup_data_ibtrack, latlon_param_TC_parmembre_pour1ech#, recup_size_ibtrack_parTC

#####################################################################
# Initializations
#####################################################################


# Récupération des variables dans la namelist
fd=open('Hurr-Track_ana.nam','r')
data_nam=np.loadtxt('Hurr-Track_ana.nam',dtype='str')[:]

annee=data_nam[0]
mois=data_nam[1]
jour=data_nam[2]
heure=data_nam[3]

print(annee, mois, jour, heure)
  
list_tc=[]
nbr_tc_tot=0

liste_bassin=["NA", "EP", "WP", "SP", "SI", "NI"]

f=open('nbcyclone.txt','w')
for ibbassin in liste_bassin :
  nbr_tc, list_tc=nombre_TC(annee, mois, jour, heure, ibbassin)
  print(nbr_tc, list_tc)
  if nbr_tc != 0 :
    f.write(str(nbr_tc)+','+ibbassin+',')
    for tc in list_tc:
     f.write(tc+',')
    f.write("\n")
  nbr_tc_tot=nbr_tc+nbr_tc_tot
if nbr_tc_tot == 0 :
  f.write(str(0))
f.close()
