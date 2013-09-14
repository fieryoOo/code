#!/bin/csh
if ( $#argv != 2 )then
echo "usage: Plot_noise_square_map.csh [Noise_square_path] [period]"
exit
endif

cd $argv[1]
rm -f 'noise_square_'$argv[2]'sec.txt'
set range=`echo $argv[2] | awk '{print 0.2*$1}'`
foreach file (`ls *_avgd_noise_square`)
set amp=`awk -v r=$range -v p=$argv[2] 'NR>1 && $1>p-r && $1<p+r' $file | awk -v p=$argv[2] 'BEGIN{sum=0;weight=0}{weight+=1/($1-p)**2;sum+=$2/($1-p)**2}END{print sum/weight}'`
set long=`head -1 $file | sed s/'\+'/''/g | awk '{print $2}'`
set long=`echo $long | awk '{if($1>=0){print $1};if($1<0){print $1+360}}'`
set lati=`head -1 $file | sed s/'\+'/''/g | awk '{print $3}'`
echo $long $lati $amp >> 'noise_square_'$argv[2]'sec.txt'
end

/home/tianye/code/Script/GMT/C_plot_travel_positive 'noise_square_'$argv[2]'sec.txt' /utera/tianye/data_ASN_TA/region_TA
/home/tianye/code/Script/GMT/TXT2CPT_vel_6decimal 'noise_square_'$argv[2]'sec.txt.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_region_contour 'noise_square_'$argv[2]'sec.txt.HD.cpt' 'noise_square_'$argv[2]'sec.txt.HD' /utera/tianye/data_ASN_TA/region_TA
