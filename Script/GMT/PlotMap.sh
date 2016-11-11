#!/bin/bash

if [ $# != 1 ] && [ $# != 2 ] && [ $# != 3 ]; then
	echo "Usage: "$0" [file-in] [plot-sta (optional, 0=no-default, 1=circle, 2=circle-with-color, 3=circle-with-name, 4=circle_only)] [fcpt (optional)]"
	exit
fi

# input
fin=$1

# correct longitudes
idel=0
fintmp=.${fin}.PlotMap.temp
ftobedeleted[idel]=$fintmp; let idel++
awk '{lon=$1; if(lon<0.){lon+=360.} print lon, $2, $3, $4}' $fin > $fintmp

# decide region
cloc=`awk 'BEGIN{lonsum=0; latsum=0}{lonsum+=$1; latsum+=$2;}END{print lonsum/NR, latsum/NR}' $fintmp`
clon=`echo $cloc | awk '{print $1}'`; clat=`echo $cloc | awk '{print $2}'`
stds=`awk -v clon=$clon -v clat=$clat 'BEGIN{lonstd=0; latstd=0}{lonstd+=($1-clon)**2; latstd+=($2-clat)**2}END{print (lonstd/(NR-1))**0.5, (latstd/(NR-1))**0.5}' $fintmp`
slon=`echo $stds | awk '{print $1}'`; slat=`echo $stds | awk '{print $2}'`
lonmin=`awk -v clon=$clon -v slon=$slon 'BEGIN{lonmin=360.}{if(lonmin>$1)lonmin=$1}END{lonmin+=(lonmin-clon)*0.05; lonm=clon-slon*2.; if(lonm>lonmin){print lonm}else{print lonmin}}' $fintmp`
lonmax=`awk -v clon=$clon -v slon=$slon 'BEGIN{lonmax=-360.}{if(lonmax<$1)lonmax=$1}END{lonmax+=(lonmax-clon)*0.05; lonm=clon+slon*2.; if(lonm<lonmax){print lonm}else{print lonmax}}' $fintmp`
#awk -v clon=$clon -v slon=$slon 'BEGIN{lonmax=-360.}{if(lonmax<$1)lonmax=$1}END{print lonmax,clon; lonmax+=(lonmax-clon)*0.05; lonm=clon+slon*2.; print lonmax, lonm}' $fintmp
latmin=`awk -v clat=$clat -v slat=$slat 'BEGIN{latmin=90.}{if(latmin>$2)latmin=$2}END{latmin+=(latmin-clat)*0.05; latm=clat-slat*2.; if(latm>latmin){print latm}else{print latmin}}' $fintmp`
latmax=`awk -v clat=$clat -v slat=$slat 'BEGIN{latmax=-90.}{if(latmax<$2)latmax=$2}END{latmax+=(latmax-clat)*0.05; latm=clat+slat*2.; if(latm<latmax){print latm}else{print latmax}}' $fintmp`
#REG=`echo $cloc $stds | awk '{print "-R"$1-$3*2."/"$1+$3*2."/"$2-$4*2."/"$2+$4*2.}'`
REG=-R${lonmin}/${lonmax}/${latmin}/${latmax}
#REG=-R1./30./0.5/6.
xmrk=`echo ${lonmin} ${lonmax} | awk '{printf "%.1f",($2-$1)/5.}'`
ymrk=`echo ${latmin} ${latmax} | awk '{printf "%.1f",($2-$1)/5.}'`
xtic=`echo $xmrk | awk '{print $1/5.}'`
ytic=`echo $ymrk | awk '{print $1/5.}'`
isGeo=false
if [ `echo $lonmin $lonmax $latmin $latmax | awk '{if($1>=-180&&$1<=360&&$2>=-180&&$2<=360&&$3>=-90&&$3<=90&&$4>=-90&&$4<=90){print 1}else{print 0}}'` == 1 ]; then
	isGeo=true
fi
echo "region = "$REG" isGeo = "$isGeo

# compute surface
res=0.1
ts=0.2
bdis=30.

fbmtmp1=${fintmp}.bm1
ftobedeleted[idel]=$fbmtmp1; let idel++
fbmtmp2=${fintmp}.bm2
ftobedeleted[idel]=$fbmtmp2; let idel++

blockmean $fintmp $REG -I$bdis'km' > $fbmtmp1
blockmean $fbmtmp1 $REG -F -I$bdis'km' > $fbmtmp2

fgrd=${fintmp}.grd
ftobedeleted[idel]=$fgrd; let idel++

surface $fbmtmp2 -T$ts -G$fgrd -I$res $REG

# plotting type
plot_type=0
if [ $# -ge 2 ]; then plot_type=$2; fi

# cpt file
if [ $# == 3 ]; then
	fcpt=$3
else
	TXT2CPT $fin
	fcpt=${fin}.cpt
	#ftobedeleted[idel]=$fcpt; let idel++
fi

# gmt
gmtset HEADER_FONT_SIZE 15
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 10
# plot starts
if [ $isGeo == ture ]; then
	SCA=-JN$clon/6i
else
	SCA=-JX16
fi
psout=${fin}.ps
pwd | psxy -H $REG $SCA -X2. -Y8. -P  -V -K  > $psout

if [ $plot_type != 4 ]; then
	sms=`echo $res | awk '{print $1/5.0}'`
	fgrds=${fintmp}.grds
	ftobedeleted[idel]=$fgrds; let idel++
	grdsample $fgrd -Q -G$fgrds $REG -I$sms
	grdimage $SCA $REG $fgrds -C$fcpt -O -K >> $psout
fi

if [ $isGeo == ture ]; then
	pscoast $SCA $REG -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 >> $psout
	dirhead=/projects/yeti4009/code/Programs/head
	if [ ! -e $dirhead ]; then dirhead=/home/tianye/code/Programs/head; fi
	psxy ${dirhead}/wus_province_II.dat $SCA $REG -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
	psxy ${dirhead}/platebound.gmt $SCA $REG -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
fi
pwd | psxy -R -J -Ba${xmrk}f${xtic}/a${ymrk}f${ytic}WeSn:."$fin": -O -K >> $psout
#pwd | psxy -R -J -Ba3f1/a3f1WeSn:."$fin": -O -K >> $psout
psscale  -C$fcpt -D6.2/-1./12.4/0.5h -O -K >> $psout

# stations
if [ ${plot_type} == 1 ]; then
	psxy $fintmp $REG $SCA -Sc.35 -W3,white -O -K >> $psout
	psxy $fintmp $REG $SCA -Sc.4 -W3,black -O -K >> $psout
elif [ ${plot_type} -gt 1 ]; then
	psxy $fintmp $REG $SCA -Sc.38 -C$fcpt -W3,black -O -K >> $psout
	#psxy $fintmp $REG $SCA -Ss.14 -C$fcpt -O -K >> $psout
fi

# texts
if [ ${plot_type} == 3 ]; then
	awk '{print $1,$2,"8 0.0 7 CB",$4}' $fintmp | pstext -R -J -O -K -Wwhite -G0 >>  $psout
fi

elon=245.123; elat=41.160
echo $elon $elat | psxy -R -J -Sc1. -Gdarkgray -W8,white -O -K >> $psout
# plot ends
pwd | psxy -R -J -O >> $psout
echo $psout

# clean up
echo ${ftobedeleted[@]} | xargs rm -f

