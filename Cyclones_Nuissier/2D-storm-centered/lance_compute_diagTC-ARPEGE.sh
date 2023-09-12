#!/usr/bin/sh

. /home/ormieresl/.epygram/profile

. /home/ormieresl/.vortexrc/profile

# ECCODES DEF
export ECCODES_SAMPLES_PATH=/usr/share/eccodes/samples
export ECCODES_DEFINITION_PATH=/usr/share/eccodes/definitions

maison='/home/ormieresl/Cyclones_Nuissier_v1/2D-storm-centered/'
tempofic='/cnrm/recyf/Work/temp/ormieresl/2D-storm_temp'
pathlola='/cnrm/recyf/Data/users/ormieresl/cyclones_nuissier/2D-storm_plots'

mkdir $tempofic
mkdir $pathlola
cd $tempofic


export PYTHONPATH=$maison/DEVPYTHONTC:$PYTHONPATH  # Outils Python développés pour caractérisation des TCs

model='arpege'
domaine='GLOB01'
#domaine='glob025' #test lola
bassin='CUSTOM-HN'
area='4dvarfr'

lvideo=false



annee=$1
mois=$2
jour=$3
heure=$4
echmax=$5
xp=$6
lrecup=$7
tcname=$8



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
  latn=35.
  lats=15.
  lono=-100.
  lone=-70.
fi



# Récupération du fichier Best-Track
#cp /cnrm/coope/nuissier/DEV-PNT-AROME-OM/BTdata/BT_HISTORY_${tcname} .
cp /cnrm/coope/nuissier/PourLola/track/BT_HISTORY_${tcname} .



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
  $maison/recup_ARPxp-Vortex.sh $annee $mois $jour $heure $model $area $domaine $xp $echmax
fi  

# Lancement du script principal

python3 $maison/Hurr-compute-diags-tc.py


if [ $lvideo = true  ]
  then 	
# Creation des repertoires de visu

#  pathsite="/home/faure/outremer/demo/photos/CASE-STUDIES/CYCLONES/${tcname}/${annee}/${annee}${mois}${jour}${heure}/arpege"
  pathsite="/home/ormieresl/Cyclones_Nuissier/2D-storm-centered/videos/${tcname}/${annee}/${annee}${mois}${jour}${heure}"

  
  mkdir -p ${pathsite}/${xp}

  x=1; for i in xp_OPER-WIND-CENTERED*.png; do for j in 0 1 2 3 4 5 6 7;do counter=$(printf %04d $x); ln -s "$i" img"$counter".png; x=$(($x+1)); done;done
  ffmpeg -framerate 20 -i img%04d.png -c:v libvpx-vp9 -pix_fmt yuva420p -lossless 1 CYCLONE-MOVIES-VORTEX-WIND.webm
  mv CYCLONE-MOVIES-VORTEX-WIND.webm ${pathsite}/${xp}

  #x=1; for i in  VIEW-CB*.png; do for j in 0 1 2 3 4 5 6 7;do counter=$(printf %04d $x); ln -s "$i" img"$counter".png; x=$(($x+1)); done;done
  #ffmpeg -framerate 20 -i img%04d.png -c:v libvpx-vp9 -pix_fmt yuva420p -lossless 1 CYCLONE-MOVIES-VORTEX-CB.webm
  #mv CYCLONE-MOVIES-VORTEX-CB.webm ${pathsite}/Storm-Centered-Plot/XY/${xp}


fi
mv $tempofic/OPER-WIND-CENTERED-*.png $pathlola/.

cd ..
rm -rf $tempofic

