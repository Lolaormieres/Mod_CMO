#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 10:32:21 2023

@author: ormieresl
"""

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

dir_csv = '/home/ormieresl/Routines/DF_BOUEES_mask/longerperiod/'
dir_im = '/cnrm/recyf/Data/users/ormieresl/plot_sst_IC/'
nb_level=12
axes = plt.gca()
liste_param=['Xp.OML','RefnoOML']
#liste_color=["black","steelblue","aquamarine"]
liste_color=["palegreen","gray"]

zone=['eurat','med','tropiques','nordat','hn20','hs20','glob']
zone_name=['Eurat','Mediteranean','Tropics', 'North Atlantic', 'North Hemisphere', 'South Hemisphere', 'Global' ]


expe=['GN84','GKCJ']
for izone in range(len(zone)):
    print(izone, zone[izone])
    df_tot = pd.read_csv(dir_csv+str(expe[0])+'_'+str(zone[izone])+'_prevision_bouees_flag_ech_maskcorr_septembre_octobre.csv')
    df_tot.head()
    df_tot.dtypes
    df_tot=df_tot[510:] #without august
    df_tot_ref = pd.read_csv(dir_csv+str(expe[1])+'_'+str(zone[izone])+'_prevision_bouees_flag_ech_maskcorr_septembre_octobre.csv')
    df_tot_ref.head()
    df_tot_ref.dtypes
  
    #df_ech.columns.values GET name of columns
    df_tot_cmo=df_tot.loc[((df_tot['expe']==expe[0]) )][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees_tot','Nb_bouees_an','diff_Nb']]
    df_tot_ref=df_tot_ref.loc[((df_tot_ref['expe']==expe[1]) )][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees_tot','Nb_bouees_an','diff_Nb']]


    df_ech_cmo=df_tot_cmo.loc[( df_tot['date_xp']=='02/09/22 00:00')][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees_tot','Nb_bouees_an','diff_Nb']]
    df_ech_ref=df_tot_ref.loc[( df_tot_ref['date_xp']=='02/09/22 00:00')][['date_valid','date_xp','expe','level','ech','biais_obs','eqm','Nb_bouees_tot','Nb_bouees_an','diff_Nb']]
    valna_cmo=df_tot_cmo
    valna_ref=df_tot_ref
    valnonanew_cmo=valna_cmo.rename(columns={"var": "niveau"})
    valnonanew_ref=valna_ref.rename(columns={"var": "niveau"}) 
    
    fig, axes = plt.subplots(nrows = 1,
                                ncols = 1, 
                                sharex = False, figsize=(12,9),dpi=150)

    sns.lineplot(data=valnonanew_cmo,x="ech",y="biais_obs",label=liste_param[0],color=liste_color[0],linewidth=2)
    sns.lineplot(data=valnonanew_ref,x="ech",y="biais_obs",label=liste_param[1],color=liste_color[1],linewidth=2)
    axes.set_xlabel("Forecast range (UTC)",fontsize=24)
    axes.tick_params(axis='x', labelsize=20)
    axes.set_ylabel("Biais (°C)",fontsize=24)
    axes.tick_params(axis='y', labelsize=20)
      
    axes.xaxis.set_ticks(range(6,102,6))

    axes.set_title(str(zone_name[izone]), fontsize = 33, y= 1.01)
  
    
    axes.grid(alpha=0.8)
    # plt.title(name[jlev])
    plt.savefig(dir_im+zone[izone]+"Plot_IC_REF_OML_eval_bouees_sept_oct.png")
plt.cla()
