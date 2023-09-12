#!/usr/bin/sh

. /home/nuissier/.epygram/profile

. /home/nuissier/.vortexrc/profile

# ECCODES DEF
export ECCODES_SAMPLES_PATH=/usr/share/eccodes/samples
export ECCODES_DEFINITION_PATH=/usr/share/eccodes/definitions

export PYTHONPATH=/cnrm/coope/nuissier/PourTom:$PYTHONPATH  # Outils Python développés pour caractérisation des TCs

model='arome'
domaine='CARAIB0025'
bassin='CUSTOM-HN'
area='antilles'
#xp='GK5W'  # Une expé rejeu avec + de params de sortie
#xp='GGIM' # 
xp='GGP0'  # Une expé rejeu pour FIONA

tcname='FIONA'
lrecup=false
lvideo=true


annee=$1
mois=$2
jour=$3
heure=$4
echmax=$5



cd tmp
rm -f RR*png VIEW-RR*png *-STORM*png Vortex*png KINEMATIC*png img*.png *.webm


if [ $bassin = 'CUSTOM-HS' ]
  then
  latn=0.
  lats=-30.
  lono=35.
  lone=90.
fi

if [ $bassin = 'CUSTOM-HN' ]
  then
  latn=40.
  lats=8.
  lone=-35.
  lono=-90.
fi



# Récupération du fichier Best-Track
cp /cnrm/coope/nuissier/DEV-PNT-AROME-OM/BTdata/BT_HISTORY_${tcname} .


# Préparation namelist pour le tracker python
cat << EOF > Hurr-Track.nam
$annee
$mois
$jour
$heure
$model
$domaine
$area
$bassin
$tcname
$echmax
EOF


if [ $lrecup = true ]
  then
  /d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/PourTom/recup_PEAROMxp-Vortex.sh $annee $mois $jour $heure $model $area $domaine $xp $echmax
fi  


# Lancement du script principal

python3 /d0/data2/faure/olivier/appli-pnt-aro-antilles/EXTRA-DOM/CH/PourTom/Hurr-compute-diags-pe.py



if [ $lvideo = true  ]
  then 	
# Creation des repertoires de visu

  pathsite="/home/faure/outremer/demo/photos/CASE-STUDIES/CYCLONES/${tcname}/${annee}/${annee}${mois}${jour}${heure}/pearo2.5km"

  
  mkdir -p ${pathsite}/Storm-Centered-Plot/XY/${xp}

  x=1; for i in  VIEW-WIND*.png; do for j in 0 1 2 3 4 5 6 7;do counter=$(printf %04d $x); ln -s "$i" img"$counter".png; x=$(($x+1)); done;done
  ffmpeg -framerate 20 -i img%04d.png -c:v libvpx-vp9 -pix_fmt yuva420p -lossless 1 CYCLONE-MOVIES-VORTEX-WIND.webm
  mv CYCLONE-MOVIES-VORTEX-WIND.webm ${pathsite}/Storm-Centered-Plot/XY/${xp}


fi



