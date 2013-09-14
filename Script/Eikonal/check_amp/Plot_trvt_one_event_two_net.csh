#!/bin/csh
if ( $#argv != 6 )then
echo "usage: Plot_trvt_two_net [sac_path_1] [sac_path_2] [period] [event] [region_infile] [v1 or v2]"
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

if ( ! -e $path1'/'$event'.ph.txt_'$qlt )then
echo "no ph "$qlt" file found for event "$event' in '$path1
exit
endif
if ( ! -e $path2'/'$event'.ph.txt_'$qlt )then
echo "no ph "$qlt" file found for event "$event' in '$path2
exit
endif

/home/tianye/code/Script/GMT/C_plot_travel $path1'/'$event'.ph.txt_'$qlt $argv[5]
/home/tianye/code/Script/GMT/C_plot_travel $path2'/'$event'.ph.txt_'$qlt $argv[5]
/home/tianye/code/Script/GMT/TXT2CPT_vel_nodecimal $path1'/'$event'.ph.txt_'$qlt'.HD'
/home/tianye/code/Script/GMT/TXT2CPT_vel_nodecimal $path2'/'$event'.ph.txt_'$qlt'.HD'

/home/tianye/code/Script/GMT/C_plot_kernel_input_two_net $cpt_path'/'$event'.ph.txt_'$qlt'.HD.cpt' $path1'/'$event'.ph.txt_'$qlt'.HD' $path2'/'$event'.ph.txt_'$qlt'.HD' $argv[5]
mv $path1'/'$event'.ph.txt_'$qlt'.HD.ps' $out_path'/'$event'.ph.'$per'sec.txt_'$qlt'.HD.ps'

