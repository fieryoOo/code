#!/bin/csh
if ($#argv != 2)then
echo "usage: Plot_hist [sac_path] [v1 or v2]"
exit
endif

cd $argv[1]
mkdir -p sta_amp_histogram
set net=`echo $argv[1] | awk -F/ '{print $(NF-1)}'`

cd '30sec_10snr_360dis/amp_diff_30s_'$argv[2]
foreach file (`ls ????_pctg sta_avg sta_std`)
cd $argv[1]
set out = $argv[1]'/sta_amp_histogram/'$file'.ps'
set sta=`echo $file | cut -d. -f1`
echo $sta
if( $file == sta_avg )then
set REG=-R-10/10/0/30
set SCL=-Jx.3c/0.2c
else if( $file == sta_std )then
set REG=-R0/20/0/30
set SCL=-Jx.3c/0.2c
else
set REG=-R-30/30/0/30
set SCL=-Jx.1c/0.2c
endif

set N = 8
set X0 = 4.1
set Y0 = `echo $N | awk '{print $1+2.5}'`
set flag = 0

foreach period (22 30 40 60 80 100)
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[2]
cp /home/tianye/code/Script/GMT/gmtdefaults4 .gmtdefaults4
gmtset HEADER_FONT_SIZE 8p
gmtset HEADER_OFFSET -0.2c
gmtset LABEL_FONT_SIZE 5p
gmtset ANNOT_FONT_SIZE_PRIMARY 5p

set temp=`awk -v sta=$sta '$1==sta {print $2"@"$3}' sta_avg_std_table`
set avg=`echo $temp | cut -d@ -f1`
set std=`echo $temp | cut -d@ -f2`
if($file == sta_avg | $file == sta_std )then
set avg=`awk 'BEGIN{a=0}{a+=$1}END{print a/NR}' $file`
endif
@ flag += 1

if ( $flag == 1)then
pshistogram $file -B6/5:."$period sec  $file  mean $avg  std $std": $REG -Z1 $SCL -W1 -L1.0p,red -F -K -Y$Y0 -X$X0 > $out
else if( $flag == 4)then
pshistogram $file -B6/5:."$period sec  $file  mean $avg  std $std": $REG -Z1 $SCL -W1 -L1.0p,red -F -K -O -Y-$N -X-`echo "2*$N" | bc` >> $out
else if( $flag == 6)then
pshistogram $file -B6/5:."$period sec  $file  mean $avg  std $std": $REG -Z1 $SCL -W1 -L1.0p,red -F -O -X$N >> $out
else
pshistogram $file -B6/5:."$period sec  $file  mean $avg  std $std": $REG -Z1 $SCL -W1 -L1.0p,red -F -K -O -X$N >> $out
endif

cd ../..
end
#exit
end
