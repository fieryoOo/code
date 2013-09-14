#!/bin/csh
if ( $#argv != 3 )then
echo "usage: Calc_azi_dec.csh [cent_sta] [sta_lst] [period]"
exit
endif

cd /home/tianye/data_center_sta/$argv[1]
set long=`grep -m1 $argv[1] $argv[2] | awk '{print $2}'`
set lati=`grep -m1 $argv[1] $argv[2] | awk '{print $3}'`
mkdir -p azi_dec_$argv[3]'_smoothed'
set file=$argv[1]'_center_normd_amp_pos_'$argv[3]'sec.txt'
cp $file azi_dec_$argv[3]'_smoothed'
cp $file'.HD' azi_dec_$argv[3]'_smoothed'
cd azi_dec_$argv[3]'_smoothed'

rm -f $argv[1]'_'$argv[3]'sec_deg_dec_amp_rsd.txt'
set deg=0
while($deg < 360)
echo 'Working on deg '$deg'...'
set deg_min=`echo $deg | awk '{print $1-10}'`
set deg_max=`echo $deg | awk '{print $1+10}'`
set deg_min_2=`echo $deg | awk '{print $1-5}'`
set deg_max_2=`echo $deg | awk '{print $1+5}'`

/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/amp_dist_angle_snr $long $lati $file $deg_min $deg_max
if( `awk '$1>600' $file'_amp_dist_deg'$deg_min'_to_deg'$deg_max | wc -l` < 3 )then
echo 'deg'$deg_min'_to_deg'$deg_max': skipped!'
@ deg += 5
continue
endif

/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/amp_dist_angle $long $lati $file'.HD' $deg_min_2 $deg_max_2

set snr=`awk '$1>600 && $1<800' $file'_amp_dist_deg'$deg_min'_to_deg'$deg_max | awk 'BEGIN{a=0}{a+=$4}END{print a/NR}'`
set dist_max=`awk 'BEGIN{a=0}{if($1>a)a=$1}END{print a}' $file'_amp_dist_deg'$deg_min'_to_deg'$deg_max`
awk -v dist=$dist_max '$1<dist {print $1,$2}' $file'.HD_amp_dist_deg'$deg_min_2'_to_deg'$deg_max_2 > temp_deg
set dec_amp_rsd=`/home/tianye/code/Programs/FIT/least_squares_line temp_deg`
echo $deg $dec_amp_rsd $snr >> $argv[1]'_'$argv[3]'sec_deg_dec_amp_rsd.txt'
@ deg += 5
end
