#!/bin/csh
if( $#argv != 9 )then
   echo "Usage: "$0" [InputMap] [MapResolution (in km)] [lon1] [lat1] [lon2] [lat2] [c1] [c2] [c3]"
   exit
endif

set ftmp = temp_PathAvg
#project to pq coordinates
set wd = `echo $argv[2] | awk '{print $1*3}'`
set dist = `project -G100000 -C$argv[3]'/'$argv[4] -E$argv[5]'/'$argv[6] -Q -S | awk 'NR==2{print $3}'`
set rend = `echo $wd $dist | awk '{print $1+$2}'`
project $argv[1] -C$argv[3]'/'$argv[4] -E$argv[5]'/'$argv[6] -L-$wd'/'$rend -Q -S -W-$wd'/'$wd -Fpqz > $ftmp
if( `more $ftmp | wc -l` < 30 ) then
   echo "No enough data points around the path!"
   exit
endif

#blockmean
set REG = '-R-'$wd'/'$rend'/-'$wd'/'$wd
set bdis = `echo $argv[2] | awk '{print $1/3.}'`
blockmean $ftmp $REG -I$bdis'km' > $ftmp'.txt1'
blockmean $ftmp'.txt1' $REG -F -I$bdis'km' > $ftmp
rm -f $ftmp'.txt1'

#surface and compute average
set ts = 0.2
set nseg = `echo $dist $argv[2] | awk '{printf "%.0f",$1*5/$2}'`
if( $nseg == 0 ) then
   echo "Increase Map resolution!"
   exit
endif
set res = `echo $dist $nseg | awk '{print $1/$2}'`
surface $ftmp -T$ts -G$ftmp'.tomo' -I$res $REG
#set wavg = `grd2xyz $ftmp'.tomo' $REG | awk -v dist=$dist -v res=$res '$1>-res/2 && $1<dist+res/2 && $2<res/2 && $2>-res/2' | awk -v c1=$argv[7] -v c2=$argv[8] -v c3=$argv[9] 'BEGIN{a=0;weit=0}{weight=1./(c1+c2*sqrt($3)+c3*$3);a+=$3*weight;weit+=weight}END{print a/weit}'`
set avg = `grd2xyz $ftmp'.tomo' $REG | awk -v dist=$dist -v res=$res '$1>-res/2 && $1<dist+res/2 && $2<res/2 && $2>-res/2' | awk 'BEGIN{a=0}{a+=$3}END{print a/NR}'`

#clean up and output
rm -f $ftmp $ftmp'.tomo'
echo $avg
