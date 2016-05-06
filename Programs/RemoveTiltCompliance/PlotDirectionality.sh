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
echo $_llon $_ulat $_title
   echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
   #echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "12. 0. 20 LT", $0}' | pstext -R -J -O -K -N >> $psout
   let ifile++
}

PlotSingleDiagram() {
	local _fin=$1
	local _t=$2
	local _psout=$3
	local _fdata=${_fin}_gmt_data_${_t}
	awk -v t=$_t '$NF==t' $_fin | awk '{azi=$1; for(i=2; i<=NF; i++){print azi,i+1,$i}}' > $_fdata
	if [ `more $_fdata | wc -l` == 0 ]; then return; fi
	#TXT2CPT $_fdata
	#local _fcpt=${_fdata}.cpt
	local _fcpt=Coherence.cpt
	local _fgrd=directionality_${_t}.grd

	#rm -f $_psout
	#pwd | psxy -Rg -JX1 -P -K > $_psout
	local _REG=-R0/360/2.5/9.5
	local _SCA=-JPa${ddiag}c
	#xyz2grd -R $fdata -G${_fgrd} -I1/0.1
	surface $_REG $_fdata -G${_fgrd} -I0.5/0.05 -T0.2
	grdimage -R $_SCA ${_fgrd} -C$_fcpt -O -K >> $_psout
	#psxy $fdata -R -J -Sc.15 -W1 -C$fcpt -A -O -K >> $_psout
	psbasemap $_REG $_SCA -Bg60/g2:"aa".: -P -O -K >> $_psout
	#pwd | psxy -R -J -O >> $_psout

	#rm	-f $_fdata $_fcpt $_fgrd
}


### main ###
stalst=/work2/tianye/ASN_OBS_Denoise/stations/station_all.lst
lon0=-131; lat0=39
REG=-R${lon0}/-120/${lat0}/50
#SCA=-JN-125.5/6i
# kmperdeg_lon = 78.8463, kmperdeg_lat = 111.1415
rdiag=1.5; ddiag=`echo $rdiag| awk '{print $1*2.}'`
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

date=2011.DEC.25
# produce a frame for each time window
fin0=ft_${date}.I05D.BHZ.SAC_RDirect
for t in `awk '{print $NF}' $fin0 | uniq`; do
	psout=${pdir}/${date}_${t}.ps; cp $psbase $psout
	for sta in J23A G03A J47A M07A J49A J65A I03D I05D; do
		fin=ft_${date}.${sta}.BHZ.SAC_RDirect
		loc=`awk -v sta=$sta '$1==sta{print $2,$3}' $stalst | head -n1`
		Slon=`echo $loc | awk -v lon0=$lon0 -v cperd=${cperd_lon} -v rdiag=$rdiag '{dlon=$1-lon0; if(dlon>180){dlon-=360.} print dlon*cperd-rdiag}'`
		Slat=`echo $loc | awk -v lat0=$lat0 -v cperd=${cperd_lat} -v rdiag=$rdiag '{print ($2-lat0)*cperd-rdiag}'`
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

convert -quality 100 -delay 300 -loop 0 -scale 100% `ls ${pdir}/${date}_*.ps | awk -F'_' '{print $0,$NF}' | sort -g -k2 | awk '{print $1}'` ${fin}.gif

echo ${fin}.gif
\rm -r $pdir
