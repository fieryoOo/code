#!/bin/csh
if ($#argv != 4)then
echo "usage: Plot_amp_std [sac_path] [sta_loc_infile] [back_g_event] [v1 or v2]"
exit
endif

cd $argv[1]
#set table_file='/home/tianye/data_Eikonal/SAC_XR/amp_diff_avg_table_XR'
set net=`echo $argv[1] | awk -F/ '{print $NF}'`
if ( ! `echo $net | awk '{print length($1)}'` )set net=`echo $argv[1] | awk -F/ '{print $(NF-1)}'`
foreach period (30 50 80)
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[4]
rm -f sta_loc_amp_std_color_$period's'
foreach sta_info (`awk '{print $1,$5}' sta_avg_std_table | sed s/'_pctg'/''/g | sed s/' '/'@'/g`)
set sta=`echo $sta_info | cut -d@ -f1`
set amp_std=`echo $sta_info | cut -d@ -f2`

set temp=`grep $sta $argv[2] | sed s/'\+'/''/g | awk '{a=$2;if(a<0){a+=360};print a,$3}'`
set long = `echo $temp | awk '{print $1}'`
set lati = `echo $temp | awk '{print $2}'`
if(`awk -v long=$long -v lati=$lati '(long-$1)^2+(lati-$2)^2<0.00001' ../$argv[3]'_am.txt_'$argv[4] | wc -l` == 0)continue
echo $temp $amp_std >> sta_loc_amp_std_color_$period's'
end
/home/tianye/code/Script/GMT/C_plot_travel ../$argv[3]'_am.txt_'$argv[4] /home/tianye/data_Eikonal/SAC_XR/region_XR
csh /home/tianye/code/Script/GMT/TXT2CPT_vel ../$argv[3]'_am.txt_'$argv[4]'.HD'
set temp_scl=`grep $period /home/tianye/data_Eikonal/SAC_TA/amp_std_scl.txt | awk '{print $2"@"$3}'`
#echo $period $temp_scl
set scl_b=`echo $temp_scl | cut -d@ -f1`
set scl_m=`echo $temp_scl | cut -d@ -f2`
set scl_e=`awk 'BEGIN{a=0}{if(a<$3){a=$3}}END{printf "%.0f",a}' sta_loc_amp_std_color_$period's'`

echo "0 255 255 255 "$scl_b" 255 255 255\
"$scl_b" 255 255 255 "$scl_m" 255 0 0\
"$scl_m" 0 0 0 "$scl_e" 0 0 0\
B 255 255 255\
F 0 0 0\
N 128 128 128\" > sta_loc_amp_std_color_$period's.cpt'

#set gray_red_cpt=`ls /home/tianye/data_Eikonal/SAC_XR/$period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[4]'/gray_red_'*'steps.cpt'`
/home/tianye/code/Script/GMT/C_plot_amp_diff ../$argv[3]'_am.txt_'$argv[4]'.HD.cpt' ../$argv[3]'_am.txt_'$argv[4]'.HD' sta_loc_amp_std_color_$period's.cpt' sta_loc_amp_std_color_$period's' /home/tianye/data_Eikonal/SAC_XR/region_XR
mkdir -p /home/tianye/data_Eikonal/results/raw_am_plus_sta_shift/$argv[3]
mv sta_loc_amp_std_color_$period's.ps' /home/tianye/data_Eikonal/results/raw_am_plus_sta_shift/$argv[3]/$argv[3]'_'$net'_'$period's_std_'$argv[4]'.ps'
cd ../..
#exit
end
