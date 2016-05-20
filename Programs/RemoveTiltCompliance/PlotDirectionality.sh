#!/bin/bash

PlotText() {
   local _REG=$1
   local _text=$2
   local _XS=$3
   local _YS=$4

   #local _lb=('a' 'b' 'c' 'd' 'e' 'f')
   #local _title=`echo ${_text} | awk -v lb=${_lb[ifile]} '{print "("lb")  "$0}'`
   local _title=$_text
   echo ${_title}
   local _llon=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v xs=$_XS '{print $1+xs*($2-$1)}'`
   local _ulat=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v ys=$_YS '{print $4+ys*($4-$3)}'`
   echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
   #echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "12. 0. 20 LT", $0}' | pstext -R -J -O -K -N >> $psout
   let ifile++
}

PlotSingleDiagram() {
	local _fin=$1
	local _t=$2
	local _psout=$3
	local _fdata=${_fin}_gmt_data_${_t}
	#awk -v t=$_t '$NF==t' $_fin | awk '{azi=$1; for(i=2; i<=NF; i++){print azi,i+1,$i}}' > $_fdata
	awk -v t=$_t '0.5*($1+$2)==t' $_fin | awk '{azi=$3; for(i=4; i<=NF; i++){print azi,i-2.5,$i}}' > $_fdata
	if [ `more $_fdata | wc -l` == 0 ]; then
		#echo "deleting "$_fdata
		rm -f $_fdata; return
	fi
	#TXT2CPT $_fdata
	#local _fcpt=${_fdata}.cpt
	local _fcpt=topo/Coherence.cpt
	local _fgrd=directionality_${_t}.grd

	#rm -f $_psout
	#pwd | psxy -Rg -JX1 -P -K > $_psout
	local _REG=-R0/360/1./10.
	local _SCA=-JPa${ddiag}c
	#xyz2grd -R $fdata -G${_fgrd} -I1/0.1
	surface $_REG $_fdata -G${_fgrd} -I0.5/0.05 -T0.2
	grdimage -R $_SCA ${_fgrd} -C$_fcpt -O -K >> $_psout
	#psxy $fdata -R -J -Sc.15 -W1 -C$fcpt -A -O -K >> $_psout
	gmtset FRAME_PEN 0.5p
	psbasemap $_REG $_SCA -Bg60/g2:"aa".: -P -O -K >> $_psout
	gmtset FRAME_PEN 1.0p
	#pwd | psxy -R -J -O >> $_psout

	#rm -f $_fdata $_fcpt $_fgrd
	rm	-f $_fdata $_fgrd
}


### main ###
if [ $# != 1 ]; then
	echo "Usage: "$0" [date (2012.MAR.1)]"
	exit
fi

gmtset FRAME_PEN 1.0p
#date=2012.MAR.5
date=$1
month=$(echo $date | cut -d. -f1,2)
echo $date $month

stalst=/work2/tianye/ASN_OBS_Denoise/stations/station_all.lst
lon0=-131; lat0=39
REG=-R${lon0}/-120/${lat0}/50
#SCA=-JN-125.5/6i
# kmperdeg_lon = 78.8463, kmperdeg_lat = 111.1415
rdiag=0.6; ddiag=`echo $rdiag| awk '{print $1*2.}'`
cperd_lon=1.34; cperd_lat=1.89
SCA=-Jx${cperd_lon}cd/${cperd_lat}cd
#SCA=-Jm-125.5/1.3c
psbase=.PlotDirectionality_basemap.ps
if [ ! -e $psbase ]; then
	grdimage topo/topo2.grd -Ctopo/Relief.cpt -Itopo/topo_gradient.grd -Xc -Yc $REG $SCA -K -P > $psbase
	pscoast $SCA $REG -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 >> $psbase
fi

pdir=.PlotDirectionality_tmp
mkdir -p $pdir

dir=${month}/${date}
# generate a frame for each time window
fin0=${dir}/ft_${date}.I05D.BHZ.SAC_RDirect
for t in `awk 'NR>1{print 0.5*($1+$2)}' $fin0 | uniq`; do
#for t in `awk 'NR>1{print ($1+$2)*0.5}' $fin0 | uniq | awk '$1>21500'`; do
	psout=${pdir}/${date}_${t}.ps
	if [ -e $psout ]; then continue; fi
	cp $psbase $psout
	#for sta in J23A J06A G03A J39A J37A J55A J68A M07A J25A J49A J65A I03D I05D O02D M02C L02D K04D J04D G03D F05D D04D NLWA A04D C06D; do
		#fin=${dir}/ft_${date}.${sta}.BHZ.SAC_RDirect
	for fin in ${dir}/ft_${date}.*.BHZ.SAC_RDirect; do
		sta=`echo $fin | cut -d/ -f3 | cut -d. -f4`
		loc=`awk -v sta=$sta '$1==sta{print $2,$3}' $stalst | head -n1`
		Slon=`echo $loc | awk -v lon0=$lon0 -v cperd=${cperd_lon} -v rdiag=$rdiag '{dlon=$1-lon0; if(dlon>180){dlon-=360.} print dlon*cperd-rdiag}'`
		Slat=`echo $loc | awk -v lat0=$lat0 -v cperd=${cperd_lat} -v rdiag=$rdiag '{print ($2-lat0)*cperd-rdiag}'`
		if [ $(echo $Slon'<0.' | bc -l) == 1 ] || [ $(echo $Slat'<0.' | bc -l) == 1 ]; then continue; fi
		pwd | psxy -R -J -X${Slon}c -Y${Slat}c -O -K >> $psout
		PlotSingleDiagram $fin $t ${psout}
		pwd | psxy $REG $SCA -X-${Slon}c -Y-${Slat}c -O -K >> $psout
	done
	psbasemap $REG $SCA -Ba5f1.0g1/a5f1.0g1:."": -V -O -K >> $psout
	title=${date}_${t}
	PlotText $REG "$title" 0. 0.02
	pwd | psxy -R -J -O >> $psout
	echo ${psout}
done

convert -quality 100 -delay 300 -loop 0 -scale 100% `ls ${pdir}/${date}_*.ps | awk -F'_' '{print $0,$NF}' | sort -g -k2 | awk '{print $1}'` ./${date}_RDirect.gif

echo $pdir ${fin}.gif
#\rm -r $pdir
