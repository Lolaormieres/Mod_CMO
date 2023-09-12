#Sorry for the mess!
Repertoire pour les images: /home/ormieresl/images_cnrm_repertory/
Repertoire pour les scripts: /home/ormieresl/Routines_total/


Evol_bouees_mercator_ACCORD.py 					--> trace comparaison avec bouees et mercator sur meme figure
Cartes_moyenne_diurne_clean.py 					--> trace cartes moyennes de rechauffement diurne
Cartes_sst.py 		       					--> trace cartes de sst
Cartes_moyenne_abs_optimise.py 					--> trace cartes d'eccart absolu avec Merc entre deux experiences differentes

obs/obs_guess.py 						--> Recuperer Netcdf (avec obs values, fgdepar, andeparm lon, lat etc..) -- fichier nc mis dans repertoire respectif de .fa
-------Repertoire dans lequel il y a les fichiers csv de biais et std de SST et/ou Temp en fct echeances et les .py 
1. DF_SST_mapfactor_mask_trop_corr/
DF_boucle_sst_expe_mapfactor_mask.py  				--> calcul biais de sst fct des echeances par rapport a Merc.  
2. DF_TEMP_mapfactor_mask/
DF_boucle_ech_carole_mapfactor_mask.py 				--> calcul biais de temp sur la colonne d'eau en fct des echeances par rapport a Merc. 
3. DF_BOUEES_mask/
DF_boucle_sst_expe_mapfactor_mask.py   				--> Calul biais et ecarts types de sst par rapport aux bouees (f.extract,lon,lat)



DF_plot_sst.py 			       				--> Plot evolution temporelle des biais de SST et std -- Avec csv DF_SST_mapfactor_mask_trop_corr/*csv
DF_plot_sst_expe_optimise_sensitivitytest.py			--> Fig. 6 -- Avec csv DF_SST_mapfactor_mask_trop_corr/*csv

plot_an_fg_depar.py 		          			--> Fig. 11 plot andepar et fgdepar issu de assimilation des bouees(avant recuperer les Netcdf avec /home/ormieresl/Routines/obs/obs_guess.py, 				indiquer date et experiences + varno = 4 (bouees fixes et derivantes), varno = 1 (bouees + navires), anflag=1 signifie data de bouee valide pour l'assimilation.

Score_ech_generale_allzone_expe_mapfactor_nottitle.py   	--> Fig.8 Patates des biais de temp sur la col d'eau
Score_ech_generale_allzone_expe_mapfactor_temp_Hn_stat.py 	--> Fig. 13 Temp moyennes au cours de la prevision - expe CMO & Mercator

Score_ech_bouees_hn_hs_biais_std.py   		 	  	--> Fig. 12 Eval a partir des bouees au cours de la prevision (biais + std) HN vs HS
Score_ech_bouees_domaine_biais_std.py				--> Fig. 12 Eval a partir des bouees au cours de la prevision (biais + std) domaine par domaine


-------Interval de confiance
Mean_spread_sst_previ_carole_v1.py    				--> Biais de SST et TEMPOC avec Merc + IC
Mean_spread_sst_previ_carole.py	     				--> Biais de SST et TEMPOC avec obs bouees + IC
Mean_spread_sst_previ_carole_bouees.py 				--> Biais de SST et TEMPOC avec obs bouees + IC (+exp assim bouees)
Mean_spread_sst_previ_carole_bouees_R10x4_std.py  		--> Std de SST et TEMPOC avec obs bouees + IC


Evol_sst_onegridpoint_102h.py					--> Fig 14. Evolution SST en un pt de griller (Med et tropiques) 
Evol_sst_wind_fsol_onegridpoint_102h.py				--> Pareil avec paran atmospherique en plus (vent ou flux solaire)



------ Grille etiree ARPEGE 
sgmeandm.py 							--> Calcul moyenne par domaine en corrigeant des coeff de taille de grille


------- 
/home/ormieresl/Cyclones_carole/
file shell & py : 									Pour lire/ecrire les differents parametres du fichier fa le long de la trace du cyclone
/home/ormieresl/Routines/Cyclones_caroles_csv/Plot_param_cyclones_v1.py			--> Plot differents parametres (u10, Pmin) sous le cyclone pour differentes experiences

/home/ormieresl/Cyclones_Nuissier/track ou Hovmoller ou 2D-storm-centered		--> Plot prevision cyclones versus IBtracs

