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

#,date_xp,ech,date_mer,cterm_mer,expe,var,level,biais,eqm,zone,compar,temp_exp_mean,temp_merc_mean
dir_csv='/home/ormieresl/Routines/gmapfactor/'
dir_im='/cnrm/recyf/Data/users/ormieresl/plot_sst_IC/'
nb_level=12
axes = plt.gca()
liste_param=['temp_exp_mean','temp_merc_mean']
#liste_color=["black","steelblue","aquamarine"]
liste_color=["palegreen","gold"]
name=['SFX.TEMP_OC'+str(i) for i in range(1,nb_level+1)]


zone=['eurat','med','tropiques','nordat','hn20','hs20','glob']
zone_name=['Eurat','Mediteranean','Tropics', 'North Atlantic', 'North Hemisphere', 'South Hemisphere', 'Global' ]
# zone=['eurat']

expe=['GN84']
for izone in range(len(zone)):
    print(izone, zone[izone])
    df_tot = pd.read_csv(dir_csv+str(expe[0])+'_'+str(zone[izone])+'_ech_septembre_octobre_prev-mer'+'.csv')
    df_tot.head()
    df_tot.dtypes
    #df_ech.columns.values GET name of columns
    df_tot=df_tot.loc[((df_tot['expe']==expe[0]) )][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar','temp_exp_mean','temp_merc_mean']]

    df_ech=df_tot.loc[((df_tot['var']=='SFX.TEMP_OC1') & (df_tot['date_xp']=='02/09/22 00:00'))][['date_xp','ech','date_mer','cterm_mer','expe','var','level','biais','eqm','zone','compar','temp_exp_mean','temp_merc_mean']]

    
    valna=df_tot
    valnonanew=valna.rename(columns={"var": "niveau"})
    fig, axes = plt.subplots(nrows = 1,
                                ncols = 1, 
                                sharex = False, figsize=(12,9),dpi=150)
    for jlev in range(len(name)):
        for jpar in range(len(liste_param)):
            vallev=valnonanew.loc[valnonanew.niveau == name[jlev], :].copy()
            a=vallev[liste_param[jpar]].astype(float)
            sns.lineplot(data=vallev,x="ech",y=a,label=liste_param[jpar],color=liste_color[jpar],linewidth=2)
            axes.set_xlabel("Forecast range (UTC)",fontsize=24)
            axes.tick_params(axis='x', labelsize=20)
            axes.set_ylabel("Mean temperature (°C)",fontsize=24)
            axes.tick_params(axis='y', labelsize=20)
            # axes.legend(loc="best",fontsize=19,labels=['Xp.OML SST','Xp.Mercator.instantanee'])
              
            axes.xaxis.set_ticks(range(6,102,6))
    
            axes.set_title(str(zone_name[izone]), fontsize = 33, y= 1.01)
      
            # ax.yaxis.set_label_coords(-0.1, .5) # place of labels
            # ax.xaxis.set_label_coords(0.5, -0.15)
            # axes.tick_params(width=2, length=4)
            
            axes.grid(alpha=0.8)
            # plt.title(name[jlev])
            plt.savefig("Plot_"+zone[izone]+'_'+name[jlev]+".png")
        plt.cla()
