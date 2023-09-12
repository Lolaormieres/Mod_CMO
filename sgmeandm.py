#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:33:31 2023

@author: ormieresl
"""
import numpy as np
#----------------------------------------------------------------------------------------#
def sgmeand(dta,mpf) :
  ompf = 1./mpf   # Inverse map factor
  ompf_sum = np.sum(ompf) # sigma of 1/m (weight)
  dta_norm = np.ma.divide(dta,mpf) # Global map factor normalisation
  dta_sgmean = np.sum(dta_norm) / ompf_sum # Summation and normalisation by weight
  return dta_sgmean
#----------------------------------------------------------------------------------------#
def sgmeandm(dta,mpf,mask=None) :
  if mask is not None :
     mask_init = np.ma.getmask(dta) 
     mask_inv  = np.logical_not(mask)
     mask_loc = np.ma.mask_or(mask_inv,mask_init)  
     dta_loc = np.ma.array(dta,mask=mask_loc)
     mpf_loc = np.ma.array(mpf,mask=mask_loc)
     return sgmeand(dta_loc,mpf_loc)
  else :   
     return sgmeand(dta,mpf)
#----------------------------------------------------------------------------------------#
#########version sgmeandm
# f_diff_mean = sgmeandm(f_diff.data,map_factor,mask=mask)

# ou

# f_diff_mean = sgmeandm(f_diff.data,map_factor)



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