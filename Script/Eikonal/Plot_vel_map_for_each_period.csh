#!/bin/csh
if ( $#argv != 4 ) then
echo "usage: Plot_all_vel [sac_path] [region_infile] [cpt_sac_path] [sta_infile]"
exit
endif

cd $argv[1]
foreach period ( 22 30 40 50 60 70 80 90 100 )
@ dis = $period * 12
set cpt=$argv[3]'/'$period'sec_10snr_'$dis'dis/'$period's_laplace_iso_ani_v1.1_scale_1.iso.HD.cpt'
cd $period'sec_10snr_'$dis'dis'
/home/tianye/code/Script/GMT/C_plot_travel $period's_iso_ani_v1.1' $argv[2]
/home/tianye/code/Script/GMT/C_plot_travel $period's_laplace_iso_ani_v1.1_scale_1.iso' $argv[2]

csh /home/tianye/code/Script/GMT/TXT2CPT_vel $period's_iso_ani_v1.1.HD'
csh /home/tianye/code/Script/GMT/TXT2CPT_vel $period's_laplace_iso_ani_v1.1_scale_1.iso.HD'
csh /home/tianye/code/Script/GMT/C_plot_kernel_input_sta $cpt $period's_iso_ani_v1.1.HD' $argv[4] $argv[2]
csh /home/tianye/code/Script/GMT/C_plot_kernel_input_sta $cpt $period's_laplace_iso_ani_v1.1_scale_1.iso.HD' $argv[4] $argv[2]
cd ..
end
