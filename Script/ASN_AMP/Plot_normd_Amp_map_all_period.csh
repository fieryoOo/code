#!/bin/csh
if ( $#argv != 2 )then
echo "usage: Plot_normd_Amp_map.csh [Normd_amp_path] [day or noise]"
exit
endif

cd $argv[1]
set per_file=`ls *'_avgd_'$argv[2]'normd_amp' | head -1`
foreach period (`awk 'NR>1 {print $1}' $per_file`)
rm -f $argv[2]'normd_amp_'$period'sec.txt'
echo 'Working on '$period'sec...'
foreach file (`ls *'_avgd_'$argv[2]'normd_amp'`)
#set amp=`awk -v r=$range -v p=$argv[2] 'NR>1 && $1>p-r && $1<p+r' $file | awk -v p=$argv[2] 'BEGIN{sum=0;weight=0}{weight+=1/($1-p)**2;sum+=$2/($1-p)**2}END{print sum/weight}'`
set amp=`awk -v p=$period 'NR>1 && ($1-p)**2<0.05 {print $2}' $file`
set long=`head -1 $file | sed s/'\+'/''/g | awk '{print $2}'`
set long=`echo $long | awk '{if($1>=0){print $1};if($1<0){print $1+360}}'`
set lati=`head -1 $file | sed s/'\+'/''/g | awk '{print $3}'`
echo $long $lati $amp >> $argv[2]'normd_amp_'$period'sec.txt'
end

/home/tianye/code/Script/GMT/C_plot_travel_positive $argv[2]'normd_amp_'$period'sec.txt' /utera/tianye/data_ASN_TA/region_TA
/home/tianye/code/Script/GMT/TXT2CPT_vel_4decimal $argv[2]'normd_amp_'$period'sec.txt.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_region_contour $argv[2]'normd_amp_'$period'sec.txt.HD.cpt' $argv[2]'normd_amp_'$period'sec.txt.HD' /utera/tianye/data_ASN_TA/region_TA
end
