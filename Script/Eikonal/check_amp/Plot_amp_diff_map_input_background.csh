#!/bin/csh
if ($#argv != 5)then
echo "usage: Plot_amp_diff [sac_path] [sta_loc_infile] [back_g_event] [v1 or v2] [back_g_cpt_sac_path]"
exit
endif

cd $argv[1]
#set table_file='/home/tianye/data_Eikonal/SAC_XR/amp_diff_avg_table_XR'
set net=`echo $argv[1] | awk -F/ '{print $NF}'`
if ( ! `echo $net | awk '{print length($1)}'` )set net=`echo $argv[1] | awk -F/ '{print $(NF-1)}'`
foreach period ( 30 50 80 )
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[4]
rm -f sta_loc_amp_color_$period's'
foreach sta_info (`awk '{print $1,$4}' sta_avg_std_table | sed s/'_pctg'/''/g | sed s/' '/'@'/g`)
set sta=`echo $sta_info | cut -d@ -f1`
set amp_diff=`echo $sta_info | cut -d@ -f2`

set temp=`grep $sta $argv[2] | sed s/'\+'/''/g | awk '{a=$2;if(a<0){a+=360};print a,$3}'`
set long = `echo $temp | awk '{print $1}'`
set lati = `echo $temp | awk '{print $2}'`
if(`awk -v long=$long -v lati=$lati '(long-$1)^2+(lati-$2)^2<0.00001' ../$argv[3]'_am.txt_'$argv[4] | wc -l` == 0)continue
echo $temp $amp_diff >> sta_loc_amp_color_$period's'
end
/home/tianye/code/Script/GMT/C_plot_travel ../$argv[3]'_am.txt_'$argv[4] /home/tianye/data_Eikonal/SAC_XR/region_XR
csh /home/tianye/code/Script/GMT/TXT2CPT_vel ../$argv[3]'_am.txt_'$argv[4]'.HD'
set step=`awk 'BEGIN{a=0}{if(a<sqrt($3^2)){a=sqrt($3^2)}}END{print a}' sta_loc_amp_color_$period's'`
set step=`echo "($step+0.5)*2/2" | bc`
makecpt -T-$step/$step/1 -C/home/tianye/code/Script/GMT/cpt/gray_white_red.cpt -Z > gray_red_$step'steps.cpt'
#set gray_red_cpt=`ls /home/tianye/data_Eikonal/SAC_XR/$period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[4]'/gray_red_'*'steps.cpt'`
/home/tianye/code/Script/GMT/C_plot_amp_diff $argv[5]'/'$period'sec_10snr_'$dis'dis/'$argv[3]'_am.txt_'$argv[4]'.HD.cpt' ../$argv[3]'_am.txt_'$argv[4]'.HD' gray_red_$step'steps.cpt' sta_loc_amp_color_$period's' /home/tianye/data_Eikonal/SAC_XR/region_XR
mkdir -p /home/tianye/data_Eikonal/results/raw_am_plus_sta_shift/$argv[3]
mv sta_loc_amp_color_$period's.ps' /home/tianye/data_Eikonal/results/raw_am_plus_sta_shift/$argv[3]/$argv[3]'_'$net'_'$period's_'$argv[4]'.ps'
cd ../..
#exit
end
