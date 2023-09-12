#/bin/sh

annee=$1
mois=$2
jour=$3
heure=$4
model=$5
area=$6
geometry=$7
xp=$8
echmax=$9

echo 'Debut recuperation PE Arome OM'
echo 'Timedeb='${timedeb}
echo 'Timefin='${timefin}
echo 'Timestep='${timestep}

date



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%% Récupération de tous les fichiers %%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/home/common/sync/vortex/vortex-olive/bin/vtxget.py -v --vapp=${model} --vconf=${area} --member='rangex(0,16)' --block=forecast --kind=gridpoint --origin=fcst --date=${annee}${mois}${jour}${heure} --cutoff=production --term='rangex(0-'${echmax}'-3)' --geometry=${geometry} --experiment=${xp} --model=${model} --format=grib --local=grid.[model]-forecast.[geometry:area]+00[term]-mb[member].grib


echo 'Fin de la recuperation des gribs PE Arome OM'
date

