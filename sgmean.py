#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 13:51:35 2023

@author: ormieresl
"""
import numpy as np

#----------------------------------------------------------------------------------------#
def sgmean(myfield) :

# Compute average of a field over a reduced, streched gaussian grid
# Note 1 : np is supposed to be numpy (import numpy as np)
# Note 2 : ma (masked array) divide method is used to avoid division by 0

  map_factor = myfield.geometry.map_factor_field()  # Extraction of the map_factor field from geometry object
  map_factor_data = map_factor.getdata()
  nonmasked_value = map_factor_data.count(1)  # Number of longitude points for each (pseudo)latitude circle
  
  myfield_data = myfield.getdata()
  myfield_data_norm = np.ma.divide(myfield_data,map_factor_data) # Global map factor normalisation
  myfield_data_norm_sumlat = np.sum(myfield_data_norm,1)   # Sum over longitude for each latitude circle
  myfield_data_norm_sumlat = np.ma.divide(myfield_data_norm_sumlat,nonmasked_value)  # Normalisation by the number of longitude points for each latitude circle
  Nb_lat = len(myfield_data_norm_sumlat)  # Nb of latitude circles
  myfield_data_sgmean = np.sum(myfield_data_norm_sumlat)/Nb_lat  # Sum over latitude normalised by the number of latitude circles
  return myfield_data_sgmean
#----------------------------------------------------------------------------------------#


#########version sgmean
# f_diff_mean = f_diff.data.mean()

# Il faut d'abord récupérer le map factor :

# map_factor = (f_diff.geometry.map_factor_field()).getdata()
# f_diff_mean = sgmeand(f_diff.data,map_factor)

# Pour faire une moyenne sur une zone, à la place de ce genre de lignes :

# mask = (lon >zones[key][0]) & (lon < zones[key][1]) & (lat > zones[key][2]) & (lat < zones[key][3])
# f_diff_mean = f_diff.data[mask].mean()

# Il faut mettre  :

# mask = (lon >zones[key][0]) & (lon < zones[key][1]) & (lat > zones[key][2]) & (lat < zones[key][3])
# map_factor = (f_diff.geometry.map_factor_field()).getdata()
# f_diff_mean = sgmeand2(f_diff.data,map_factor,mask)