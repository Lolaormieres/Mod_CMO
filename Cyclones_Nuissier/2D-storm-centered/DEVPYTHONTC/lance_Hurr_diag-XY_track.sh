#!/usr/bin/sh

. /home/nuissier/.epygram/profile

. /home/nuissier/.vortexrc/profile

# ECCODES DEF
export ECCODES_SAMPLES_PATH=/usr/share/eccodes/samples
export ECCODES_DEFINITION_PATH=/usr/share/eccodes/definitions

export PYTHONPATH=/home/nuissier/PYTHONDEV/TC-TOOLS:$PYTHONPATH


annee='2022'
mois='09'
jour='15'
heure='18'
model='arome'
xp='GO8P'
domaine=$1
bassin=$2
area=$3
tcname='FIONA'
lrecup=true
lvideo=true



cd tempo
rm -f RR*png VIEW-RR*png *-STORM*png Vortex*png KINEMATIC*png img*.png *.webm


if [ $bassin = 'SOOI' ]
  then
  latn=-7.2
  lats=-16.
  #lats=-25.9
  #lono=32.7
  #lone=67.6
  lone=54.
  lono=36.

fi


if [ $bassin = 'CUSTOM-HN' ]
  then
  #latn=25.
  #lats=5.
  #lono=-65.
  #lone=-20.
  latn=40.
  lats=8.
  lone=-35.
  lono=-90.
  latn3=16.8
  lats3=15.7
  lono3=-62.
  lone3=-60.9
fi

if [ $bassin = 'CUSTOM-HS' ]
  then
  latn=0.
  lats=-30.
  lono=35.
  lone=90.
  latn3=-20.7
  lats3=-21.6
  lono3=55.1
  lone3=56.
  #latn=-15.
  #lats=-20.
  #lono=170.
  #lone=180.
fi




cp /cnrm/coope/nuissier/DEV-PNT-AROME-OM/publications/PourLACy/trace-aro/legendewindmax_bw-sooi.png /cnrm/coope/nuissier/DEV-PNT-AROME-OM/publications/PourLACy/trace-aro/tempo

#cp /cnrm/coope/nuissier/DEV-PNT-AROME-OM/BTdata/BT_HISTORY_${tcname}_FAKE .
cp /cnrm/coope/nuissier/DEV-PNT-AROME-OM/BTdata/BT_HISTORY_${tcname} .

# Préparation namelist pour le tracker python
cat << EOF > Hurr-Track.nam
$annee
$mois
$jour
$heure
$model
$domaine
$latn3
$lats3
$lono3
$lone3
$bassin
$tcname
EOF


if [ $lrecup = true ]
  then
  /d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/recup_PEAROMxp-Vortex.sh $annee $mois $jour $heure $model $area $domaine $xp 	  
  #/d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/recup_MOD-Vortex.sh $annee $mois $jour $heure $model $area $domaine $xp	  
  #/d0/data2/faure/olivier/appli-pnt-aro-antilles/MDR/CH/recup_AROxp-Track.sh $annee $mois $jour $heure $model $domaine $xp $user
fi  


#python3 /d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/Hurr-Track_PNT_aro-movie-RR-V10m.py
#python3 /d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/Hurr-Track_PNT_aro-movie-RMW-Cb-Vrad10m.py
#mogrify -crop 1000x1000 Vortex-*png

python3 /d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/Hurr-kinematic-vortex-pe-RR-map.py
#python3 /d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/Hurr-kinematic-vortex-pe-WIND-map.py
#python3 /d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/Hurr-kinematic-vortex-pe-Flux.py

if [ $lvideo = true  ]
  then 	
# Creation des sous repertoire de visu

  if [ $mois = '06'  ]
    then
    ssdir="juin"
    elif [ $mois = '07' ]
    then
    ssdir="juil"
    elif [ $mois = '08' ]
    then
    ssdir="aout"
    elif [ $mois = '09' ]
    then
    ssdir="sept"
    elif [ $mois = '10' ]
    then
    ssdir="oct"
    elif [ $mois = '11' ]
    then
    ssdir="nov"
    elif [ $mois = '12' ]
    then
    ssdir="dec"
  fi

  if [ $model = 'arome' ]
     then	  
     #pathsite="/home/faure/outremer/demo/photos/CASE-STUDIES/CYCLONES/${tcname}/${annee}/${annee}${mois}${jour}${heure}/aro2.5km"
     pathsite="/home/faure/outremer/demo/photos/CASE-STUDIES/CYCLONES/${tcname}/${annee}/${annee}${mois}${jour}${heure}/pearo2.5km"
     elif [ $model = 'arpege' ]
     then	    
     pathsite="/home/faure/outremer/demo/photos/CASE-STUDIES/CYCLONES/${tcname}/${annee}/${annee}${mois}${jour}${heure}/arpege"	    
  fi    

  #mkdir -p ${pathsite}/Storm-Centered-Plot/${xp}
  mkdir -p ${pathsite}/Storm-Centered-Plot/XY/${xp}/QUANT


  x=1; for i in  VIEW-RR*.png; do for j in 0 1 2 3 4 5 6 7;do counter=$(printf %04d $x); ln -s "$i" img"$counter".png; x=$(($x+1)); done;done
  ffmpeg -framerate 20 -i img%04d.png -c:v libvpx-vp9 -pix_fmt yuva420p -lossless 1 CYCLONE-MOVIES-VORTEX-RR.webm
  mv CYCLONE-MOVIES-VORTEX-RR.webm ${pathsite}/Storm-Centered-Plot/XY/${xp}/QUANT

  #x=1; for i in  VIEW-WIND*.png; do for j in 0 1 2 3 4 5 6 7;do counter=$(printf %04d $x); ln -s "$i" img"$counter".png; x=$(($x+1)); done;done
  #ffmpeg -framerate 20 -i img%04d.png -c:v libvpx-vp9 -pix_fmt yuva420p -lossless 1 CYCLONE-MOVIES-VORTEX-WIND.webm
  #mv CYCLONE-MOVIES-VORTEX-WIND.webm ${pathsite}/Storm-Centered-Plot/XY/${xp}/QUANT


  #x=1; for i in  LATENT*.png; do for j in 0 1 2 3 4 5 6 7;do counter=$(printf %04d $x); ln -s "$i" img"$counter".png; x=$(($x+1)); done;done
  #ffmpeg -framerate 20 -i img%04d.png -c:v libvpx-vp9 -pix_fmt yuva420p -lossless 1 CYCLONE-MOVIES-VORTEX-LATENT.webm
  #mv CYCLONE-MOVIES-VORTEX-LATENT.webm ${pathsite}/Storm-Centered-Plot/XY/${xp}


fi



