

. /home/ormieresl/.epygram/profile

date

model='arpege'
domaine=$1 #GLOB025 
bassin=$2 #CUSTOM-HN or CUSTOM-HS
area=$3 #4dvarfr
echmax=$4
xp=$5
lrecup=false
ltracker=true
ldepot=false
tcname='FIONA' #IAN, FIONA does not work


annee="2022"
mois="09"
jour="21" #IAN 27
heure='00'

maison='/home/ormieresl/Cyclones_Nuissier_v1/track/'
tempofic='/cnrm/recyf/Work/temp/ormieresl/tempoarp'
pathlola='/cnrm/recyf/Data/users/ormieresl/cyclones_nuissier/track/'

mkdir $tempofic
mkdir $pathlola
cp $maison/legendewindmax*png $tempofic/.
cd $tempofic

rm -fr Vortex*png q*png prob*png track-*png *~ *~~ *~~~ img*



if [ $bassin = 'CUSTOM-HN' ] 
  then
  #latn1=35.
  #lats1=15.
  #lono1=-100.
  #lone1=-70.
  latn1=54.7 #22.0
  lats1=14.7 # instead of 9.7
  lono1=-75.3
  lone1=-51.7
fi 


if [ $bassin = 'CUSTOM-HS' ]
  then
  latn1=-10.
  lats1=-25.
  lono1=45.
  lone1=70.
  #latn1=-15.
  #lats1=-20.
  #lono1=170.
  #lone1=180.
fi


# Preparation namelist pour le tracker Aro python
cat << EOF > Hurr-Track-Aro.nam
$annee
$mois
$jour
$heure
$model
$domaine
$latn1
$lats1
$lono1
$lone1
$bassin
$tcname
$echmax
$xp
EOF


#cp /cnrm/coope/nuissier/DEV-PNT-AROME-OM/BTdata/BT_HISTORY_${tcname} .
cp /cnrm/coope/nuissier/PourLola/track/BT_HISTORY_${tcname} .


date
# Recuperation des gribs PE Arome-OM

if [ $lrecup = true ]
   then
   $maison/recup_ARPxp-Vortex.sh $annee $mois $jour $heure $model $area $domaine $xp $echmax 
fi   
 
date

# Lancement du tracker cyclones dans Arome

if [ $ltracker = true ] 
  then	
  python3 $maison/Hurr-Track_PNT_arp_gf.py
  mogrify -crop 1500x1000 ${xp}_track*png
  mv ${xp}_track*png ${pathlola}
fi  
rm -rf $tempofic

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

  mkdir -p ${pathsite}/Track

  mv ${xp}_track*png ${pathsite}/Track

fi

date

