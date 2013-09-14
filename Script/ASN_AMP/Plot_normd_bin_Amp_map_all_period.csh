#!/bin/csh
if ( $#argv != 2 )then
echo "usage: Plot_normd_bin_Amp_map.csh [Normd_bin_amp_path] [day or noise]"
exit
endif

cd $argv[1]
#foreach period (7.378866)
foreach period (`more per_file`)
/home/tianye/code/Script/GMT/C_plot_travel_positive $argv[2]'normd_'$period'.txt' /home/tianye/data_days/region_TA
/home/tianye/code/Script/GMT/TXT2CPT_vel_6decimal $argv[2]'normd_'$period'.txt.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_region_contour $argv[2]'normd_'$period'.txt.HD.cpt' $argv[2]'normd_'$period'.txt.HD' /home/tianye/data_days/region_TA
end
