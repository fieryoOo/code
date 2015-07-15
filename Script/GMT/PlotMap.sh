#!/bin/bash

if [ $# != 1 ] && [ $# != 2 ] && [ $# != 3 ]; then
	echo "Usage: "$0" [file-in] [plot-text (optional, 0=no-default, 1=yes)] [fcpt (optional)]"
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
#slon=`echo $stds | awk '{print $1}'`; slat=`echo $stds | awk '{print $2}'`
REG=`echo $cloc $stds | awk '{print "-R"$1-$3*2."/"$1+$3*2."/"$2-$4*2."/"$2+$4*2.}'`
echo "region = "$REG

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

# cpt file
if [ $# == 3 ]; then
	fcpt=$3
else
	TXT2CPT $fin
	fcpt=${fin}.cpt
	#ftobedeleted[idel]=$fcpt; let idel++
fi

# plot starts
SCA=-JN$clon/5i
psout=${fin}.ps
pwd | psxy -H $REG $SCA -X4. -Y8. -P  -V -K  > $psout

sms=`echo $ts | awk '{print $1/5.0}'`
fgrds=${fintmp}.grds
ftobedeleted[idel]=$fgrds; let idel++
grdsample $fgrd -Q -G$fgrds $REG -I$sms
grdimage $SCA $REG $fgrds -C$fcpt -Ba3f2/a3f2WeSn:."$fin": -O -K >> $psout

psxy $fintmp $REG $SCA -Sc.2 -W3,white -O -K >> $psout
psxy $fintmp $REG $SCA -Sc.25 -W3,black -O -K >> $psout

pscoast $SCA $REG -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 >> $psout

dirhead=/projects/yeti4009/code/Programs/head
if [ ! -e $dirhead ]; then dirhead=/home/tianye/code/Programs/head; fi
psxy ${dirhead}/wus_province_II.dat $SCA $REG -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
psxy ${dirhead}/platebound.gmt $SCA $REG -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
psscale  -C$fcpt -D6.2/-1./12.4/0.5h -O -K >> $psout

# texts
if [ $# -ge 2 ] && [ $2 == 1 ]; then
	awk '{print $1,$2,"8 0.0 7 CB",$4}' $fintmp | pstext -R -J -O -K -G0 >>  $psout
fi

# plot ends
pwd | psxy -R -J -O >> $psout
echo $psout

# clean up
echo ${ftobedeleted[@]} | xargs rm -f

