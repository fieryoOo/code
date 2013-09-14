#!/bin/csh
if ($#argv != 3)then
echo "usage: Statistic_avg [sac_path] [sta_loc_infile] [v1 or v2]"
exit
endif

cd $argv[1]
#set table_file='/home/tianye/data_Eikonal/SAC_XR/amp_diff_avg_table_XR'
foreach period (22 30 40 50 60 70 80 90 100)
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[3]
rm -f sta_loc_amp_color_$period's'
foreach sta_info (`awk '{print $1,$4}' sta_avg_std_table | sed s/'_pctg'/''/g | sed s/' '/'@'/g`)
set sta=`echo $sta_info | cut -d@ -f1`
set amp_diff=`echo $sta_info | cut -d@ -f2`

set temp=`grep $sta $argv[2] | sed s/'\+'/''/g | awk '{a=$2;if(a<0){a+=360};print a,$3}'`
echo $temp $amp_diff >> sta_loc_amp_color_$period's'
end
/home/tianye/code/Script/GMT/C_plot_travel /home/tianye/data_Eikonal/SAC_TA/$period'sec_10snr_'$dis'dis/'$period's_laplace_iso_ani_v1.1_scale_1.iso' /home/tianye/data_Eikonal/SAC_TA/region_TA
#csh /home/tianye/code/Script/GMT/TXT2CPT_vel /home/tianye/data_Eikonal/SAC_TA/$period'sec_10snr_'$dis'dis/'$period's_laplace_iso_ani_v1.1_scale_1.iso.HD'
#csh /home/tianye/code/Script/GMT/TXT2CPT_gray_red sta_loc_amp_color_$period's'
set step=`awk 'BEGIN{a=0}{if(a<sqrt($3^2)){a=sqrt($3^2)}}END{print a}' sta_loc_amp_color_$period's'`
set step=`echo "($step+0.5)*2/2" | bc`
makecpt -T-$step/$step/1 -C/home/tianye/code/Script/GMT/cpt/gray_red.cpt -Z > gray_red_$step'steps.cpt'
set gray_red_cpt=`ls $argv[1]'/'$period'sec_10snr_'$dis'dis/amp_diff_'$period's'$argv[3]'/gray_red_'*'steps.cpt'`
/home/tianye/code/Script/GMT/C_plot_amp_diff /home/tianye/data_Eikonal/SAC_TA/$period'sec_10snr_'$dis'dis/'$period's_laplace_iso_ani_v1.1_scale_1.iso.HD.cpt' /home/tianye/data_Eikonal/SAC_TA/$period'sec_10snr_'$dis'dis/'$period's_laplace_iso_ani_v1.1_scale_1.iso.HD' /home/tianye/data_Eikonal/station_TA.loc sta_loc_amp_color_$period's' /home/tianye/data_Eikonal/SAC_TA/region_TA $gray_red_cpt
cd ../..
#exit
end
