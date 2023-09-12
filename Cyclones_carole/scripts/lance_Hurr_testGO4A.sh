#!/bin/sh

Maison="/home/ormieresl/Cyclones_carole/"

#for jour in ${num_jour[@]}
#  do
#  for mb in `seq -f %03g 0 0`
#    do
#     member=mb$mb
#     echo $member
#  done
#done

#---------------------------
annee='2022'
heure='00'
ech_min=$1
ech_max=$2
pdt="6"
model='arpege'
domaine='glob01' #glob025
xp='GP84'
userarchi='m/mxpt/mxpt001'
#userarchi='labadie'
#userarchi='cebron'

echo $userarchi

mkdir ${Maison}Exp/${xp}
mkdir ${Maison}Exp/${xp}/${domaine}
workdir=${Maison}Exp/${xp}/${domaine}

mkdir /cnrm/recyf/Work/PEARP/Cyclones/$xp/
mkdir /cnrm/recyf/Work/PEARP/Cyclones/$xp/$domaine
patharchive='/cnrm/recyf/Work/PEARP/Cyclones/'$xp'/'$domaine'/'
cache='/cnrm/recyf/Work/temp/ormieresl/NO_SAVE/cache/vortex/arpege/4dvarfr'

cd $workdir
cp ${Maison}scripts/Hurr-Track_PNT_sst_flux.py .
cp ${Maison}scripts/fonctions2.py .
#cp ${Maison}scripts/fonctions2.pyc .
cp ${Maison}scripts/cmpt_nbcycl.py .

for mois in '09'
do
for jour in  `seq -f %02g 27 27` 
  do 
  cd $workdir
  rm -fr $workdir/tempo
  mkdir -p $workdir/tempo
  cp ${Maison}ibtracs.last3years.list.v04r00.csv $workdir/tempo
  cd $workdir/tempo
    cat << EOF > Hurr-Track_ana.nam
    $annee
    $mois
    $jour
    $heure
EOF
  python3 $workdir/cmpt_nbcycl.py
  testnbcycl=`cat nbcyclone.txt | cut -c1`
  testnbcycltout=`cat nbcyclone.txt`
  if [[ $testnbcycl = 0 ]]
  then
	  echo "youpi !" $testnbcycl " " $testnbcycltout
	  rm -f nbcyclone.txt
	  continue
  else
	  echo "chargement immediat " $testnbcycl " " $testnbcycltout
  fi

  for mb in `seq -f %03g 0 0`
   do
    member=mb$mb
# Préparation namelist pour le tracker python
    cat << EOF > Hurr-Track.nam
    $annee
    $mois
    $jour
    $heure
    $pdt
    $model
    $member
    $domaine
    $ech_min
    $ech_max
    $xp
EOF
#    python3 $workdir/Hurr-Track_PNT_sst_flux.py
    python3  $workdir/Hurr-Track_PNT_sst_flux.py
#    rm *.grib
#    rm $patharchive/*.grib
  done
  rm -f nbcyclone.txt
  rm ibtracs*
  mv *.csv $workdir
  pwd
done
done
rm -fr $workdir/tempoA
rm -rf $cache/$xp/$annee$mois$jour'T'$heure'00P/forecast/*'

exit
stop

for mois in '09' 
do
for jour in `seq -f %02g 29 30`
  do 
  cd $workdir
  rm -fr $workdir/tempo
  mkdir -p $workdir/tempo
  cp ${Maison}ibtracs.last3years.list.v04r00.csv $workdir/tempo
  cd $workdir/tempo
    cat << EOF > Hurr-Track_ana.nam
    $annee
    $mois
    $jour
    $heure
EOF
  python3 $workdir/cmpt_nbcycl.py
  testnbcycl=`cat nbcyclone.txt | cut -c1`
  testnbcycltout=`cat nbcyclone.txt`
  if [[ $testnbcycl = 0 ]]
  then
	  echo "youpi !" $testnbcycl " " $testnbcycltout
	  rm -f nbcyclone.txt
	  continue
  else
	  echo "chargement immediat " $testnbcycl " " $testnbcycltout
  fi

  for mb in `seq -f %03g 0 34`
   do
    member=mb$mb
    cd $patharchive
    sh $workdir/recup_PEArp-Track.sh $annee $mois $jour $heure $pdt $model $domaine $xp $member $ech_min $ech_max $userarchi
    cd $workdir/tempo
    ln -s $patharchive/grid.$model-forecast.$domaine+*:00.grib .
# Préparation namelist pour le tracker python
    cat << EOF > Hurr-Track.nam
    $annee
    $mois
    $jour
    $heure
    $pdt
    $model
    $member
    $domaine
    $ech_min
    $ech_max
EOF
    python3 $workdir/Hurr-Track_PNT.py
    rm *.grib
    rm $patharchive/*.grib
  done
  rm -f nbcyclone.txt
  rm ibtracs*
  mv *.csv $workdir
  pwd
done
done


for mois in '10' 
do
for jour in `seq -f %02g 01 05` `seq -f %02g 16 18` `seq -f %02g 21 31`
  do 
  cd $workdir
  rm -fr $workdir/tempo
  mkdir -p $workdir/tempo
  cp ${Maison}ibtracs.last3years.list.v04r00.csv $workdir/tempo
  cd $workdir/tempo
    cat << EOF > Hurr-Track_ana.nam
    $annee
    $mois
    $jour
    $heure
EOF
  python3 $workdir/cmpt_nbcycl.py
  testnbcycl=`cat nbcyclone.txt | cut -c1`
  testnbcycltout=`cat nbcyclone.txt`
  if [[ $testnbcycl = 0 ]]
  then
	  echo "youpi !" $testnbcycl " " $testnbcycltout
	  rm -f nbcyclone.txt
	  continue
  else
	  echo "chargement immediat " $testnbcycl " " $testnbcycltout
  fi

  for mb in `seq -f %03g 0 34`
   do
    member=mb$mb
    cd $patharchive
    sh $workdir/recup_PEArp-Track.sh $annee $mois $jour $heure $pdt $model $domaine $xp $member $ech_min $ech_max $userarchi
    cd $workdir/tempo
    ln -s $patharchive/grid.$model-forecast.$domaine+*:00.grib .
# Préparation namelist pour le tracker python
    cat << EOF > Hurr-Track.nam
    $annee
    $mois
    $jour
    $heure
    $pdt
    $model
    $member
    $domaine
    $ech_min
    $ech_max
EOF
  #  python3 $workdir/Hurr-Track_PNT.py
#    python3 $workdir/Hurr-Track_PNT_sst_flux.py
    python3 $workdir/Hurr-Track_PNT_sstsfx_flux.py

    rm *.grib
    rm $patharchive/*.grib
  done
  rm -f nbcyclone.txt
  rm ibtracs*
  mv *.csv $workdir
  pwd
done
done
