#!/usr/bin/sh


date

annee='2022'
mois='09'
jour='22'#IAN 27
heure='00'
model='arpege'
domaine=$1
bassin=$2
area=$3
echmax=$4
xp=$5
tcname='FIONA' # IAN
lrecup=true
ltracker=true
ldepot=false

maison='/home/ormieresl/Cyclones_Nuissier/Hovmoller'
cd tempoarp

cp /cnrm/coope/nuissier/DEV-PNT-AROME-OM/BTdata/BT_HISTORY_${tcname} .



# Préparation namelist pour le tracker python
cat << EOF > Hurr-Track.nam
$annee
$mois
$jour
$heure
$model
$domaine
$bassin
$xp
$tcname
$echmax
EOF


# Recuperation des gribs PE Arome-OM

if [ $lrecup = true ]
   then
   $maison/recup_ARPxp-Vortex.sh $annee $mois $jour $heure $model $area $domaine $xp $echmax
fi

date




python3 $maison/Hurr_diag-Hovmo_radius.py



if [ $ldepot = true  ]
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

  pathsite="/home/faure/outremer/demo/photos/CASE-STUDIES/CYCLONES/${tcname}/${annee}/${annee}${mois}${jour}${heure}/arpege"

  mkdir -p ${pathsite}/Hovmoller-Plots

  mv Hovmo*png ${pathsite}/Hovmoller-Plots

fi


