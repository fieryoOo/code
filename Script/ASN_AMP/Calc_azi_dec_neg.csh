#!/bin/csh
if ( $#argv != 3 )then
echo "usage: Calc_azi_dec.csh [cent_sta] [sta_lst] [period]"
exit
endif

cd /home/tianye/data_center_sta/$argv[1]
set long=`grep -m1 $argv[1] $argv[2] | awk '{print $2}'`
set lati=`grep -m1 $argv[1] $argv[2] | awk '{print $3}'`

mkdir -p azi_dec_neg_$argv[3]
set file=$argv[1]'_center_normd_amp_neg_'$argv[3]'sec.txt'
cp $file azi_dec_neg_$argv[3]
cd azi_dec_neg_$argv[3]

rm -f $argv[1]'_'$argv[3]'sec_deg_dec_amp_rsd.txt'
set deg=0
while($deg < 360)
echo 'Working on deg '$deg'...'
set deg_min=`echo $deg | awk '{print $1-10}'`
set deg_max=`echo $deg | awk '{print $1+10}'`
/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/amp_dist_angle_snr $long $lati $file $deg_min $deg_max
if( `awk '$1>500' $file'_amp_dist_deg'$deg_min'_to_deg'$deg_max | wc -l` < 10 )then
echo 'deg'$deg_min'_to_deg'$deg_max': skipped!'
@ deg += 5
continue
endif
set snr=`awk '$1>500 && $1<700' $file'_amp_dist_deg'$deg_min'_to_deg'$deg_max | awk 'BEGIN{a=0}{a+=$4}END{print a/NR}'`
awk '{print $1,$2}' $file'_amp_dist_deg'$deg_min'_to_deg'$deg_max > temp_deg
set dec_amp_rsd=`/home/tianye/code/Programs/FIT/least_squares_line temp_deg`
echo $deg $dec_amp_rsd $snr >> $argv[1]'_'$argv[3]'sec_deg_dec_amp_rsd.txt'
@ deg += 5
end
