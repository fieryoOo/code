#!/bin/bash

if [ $# != 1 ] && [ $# != 2 ] && [ $# != 3 ]; then
	echo "Usage: "$0" [file-in] [plot-sta (optional, 0=no-default, 1=circle, 2=circle-with-color, 3=circle-with-name, 4=circle_only)] [fcpt (optional)]"
	exit
fi

# input
fin=$1

idel=0
fintmp=.${fin}.PlotMap.temp
ftobedeleted[idel]=$fintmp; let idel++
cp $fin $fintmp

###
ranges=`minmax -C $fin`
lonmin=`echo $ranges | awk '{print $1}'`
lonmax=`echo $ranges | awk '{print $2}'`
latmin=`echo $ranges | awk '{print $3}'`
latmax=`echo $ranges | awk '{print $4}'`
REG=-R${lonmin}/${lonmax}/${latmin}/${latmax}
#REG=-R1./30./0.5/6.
xmrk=`echo ${lonmin} ${lonmax} | awk '{printf "%.1f",($2-$1)/5.}'`
ymrk=`echo ${latmin} ${latmax} | awk '{printf "%.1f",($2-$1)/5.}'`
xtic=`echo $xmrk | awk '{print $1/5.}'`
ytic=`echo $ymrk | awk '{print $1/5.}'`
isGeo=false
echo "region = "$REG" isGeo = "$isGeo

# compute surface
res=0.1
ts=0.2

fgrd=${fintmp}.grd
ftobedeleted[idel]=$fgrd; let idel++

surface ${fintmp} -T$ts -G$fgrd -I$res $REG

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
	#psxy $fintmp $REG $SCA -Sc.38 -C$fcpt -W3,black -O -K >> $psout
	psxy $fintmp $REG $SCA -Ss.14 -C$fcpt -O -K >> $psout
fi

# texts
if [ ${plot_type} == 3 ]; then
	awk '{print $1,$2,"8 0.0 7 CB",$4}' $fintmp | pstext -R -J -O -K -Wwhite -G0 >>  $psout
fi

# plot ends
pwd | psxy -R -J -O >> $psout
echo $psout

# clean up
echo ${ftobedeleted[@]} | xargs rm -f

