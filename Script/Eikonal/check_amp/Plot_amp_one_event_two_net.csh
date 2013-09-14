#!/bin/csh
if ( $#argv != 6 )then
echo "usage: Plot_amp_two_net [sac_path_1] [sac_path_2] [period] [event] [region_infile] [v1 or v2]"
exit
endif

set out_path=/home/tianye/data_Eikonal/results/check_ph_amp_XR
set per = $argv[3]
set event = $argv[4]
@ dis = $per * 12
set path1=`echo $argv[1]'/'$per'sec_10snr_'$dis'dis'`
set path2=`echo $argv[2]'/'$per'sec_10snr_'$dis'dis'`
set cpt_path=$path1
set qlt=$argv[6]
set label1=`echo $argv[1] | awk -F_ '{print $NF}'`
set label2=`echo $argv[2] | awk -F_ '{print $NF}'`

if ( ! -e $path1'/'$event'_am.txt_'$qlt )then
echo "no am "$qlt" file found for event "$event' in '$path1
exit
endif
if ( ! -e $path2'/'$event'_am.txt_'$qlt )then
echo "no am "$qlt" file found for event "$event' in '$path2
exit
endif

/home/tianye/code/Script/GMT/C_plot_travel $path1'/'$event'_am.txt_'$qlt $argv[5]
/home/tianye/code/Script/GMT/C_plot_travel $path2'/'$event'_am.txt_'$qlt $argv[5]
/home/tianye/code/Script/GMT/TXT2CPT_vel $path1'/'$event'_am.txt_'$qlt'.HD'
/home/tianye/code/Script/GMT/TXT2CPT_vel $path2'/'$event'_am.txt_'$qlt'.HD'

/home/tianye/code/Script/GMT/C_plot_kernel_input_region $cpt_path'/'$event'_am.txt_'$qlt'.HD.cpt' $path1'/'$event'_am.txt_'$qlt'.HD' $argv[5]
/home/tianye/code/Script/GMT/C_plot_kernel_input_region $cpt_path'/'$event'_am.txt_'$qlt'.HD.cpt' $path2'/'$event'_am.txt_'$qlt'.HD' $argv[5]
#/home/tianye/code/Script/GMT/C_plot_kernel_input_two_net $cpt_path'/'$event'_am.txt_'$qlt'.HD.cpt' $path1'/'$event'_am.txt_'$qlt'.HD' $path2'/'$event'_am.txt_'$qlt'.HD' $argv[5]
mv $path1'/'$event'_am.txt_'$qlt'.HD.ps' $out_path'/'$event'_am.'$per'sec.'$label1'.txt_'$qlt'.HD.ps'
mv $path2'/'$event'_am.txt_'$qlt'.HD.ps' $out_path'/'$event'_am.'$per'sec.'$label2'.txt_'$qlt'.HD.ps'

