#!/bin/csh
if ( $#argv != 4 )then
echo "usage: Plot_trvt_amp [sac_path] [period] [event] [region_infile]"
exit
endif

cd $argv[1]
set per = $argv[2]
set event = $argv[3]
@ dis = $per * 12
cd $per'sec_10snr_'$dis'dis'
set cpt_path='/home/tianye/data_Eikonal/SAC_XR/'$per'sec_10snr_'$dis'dis/'

if ( ! -e $event'.ph.txt_v2' )then
echo "no phase v2 file found for event "$event'!'
exit
endif
if ( ! -e $event'_am.txt_v2' )then
echo "no amp v2 file found for event "$event'!'
exit
endif
mkdir -p $event.graphs

/home/tianye/code/Script/GMT/C_plot_travel $event'.ph.txt_v2' $argv[4]
/home/tianye/code/Script/GMT/TXT2CPT_vel_nodecimal $event'.ph.txt_v2.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_region $cpt_path$event'.ph.txt_v2.HD.cpt' $event'.ph.txt_v2.HD' $argv[4]
mv $event'.ph.txt_v2.HD.ps' $event.graphs

#/home/tianye/code/Script/GMT/C_plot_travel $event'_am.txt_v2' $argv[4]
#/home/tianye/code/Script/GMT/TXT2CPT_vel $event'_am.txt_v2.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_region_contour $cpt_path$event'_am.txt_v2.HD.cpt' $event'_am.txt_v2.HD' $argv[4]
mv $event'_am.txt_v2.HD.ps' $event.graphs

grep -v '0 999' 'slow_azi_'$event'.txt.HD' | awk '{print $1,$2,$3**2}' > 'slow_term_'$event'.txt.HD'
/home/tianye/code/Script/GMT/TXT2CPT_vel_3decimal 'slow_term_'$event'.txt.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_sta $cpt_path'slow_term_'$event'.txt.HD.cpt' 'slow_term_'$event'.txt.HD' $event'_am.txt_v2' $argv[4]
mv 'slow_term_'$event'.txt.HD.ps' $event.graphs

grep -v '0 999' 'slow_azi_'$event'.txt.HD' | awk '{print $1,$2,1/$3}' > 'vel_apparent_'$event'.txt.HD'
/home/tianye/code/Script/GMT/TXT2CPT_vel 'vel_apparent_'$event'.txt.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_sta $cpt_path'vel_apparent_'$event'.txt.HD.cpt' 'vel_apparent_'$event'.txt.HD' $event'_am.txt_v2' $argv[4]
mv 'vel_apparent_'$event'.txt.HD.ps' $event.graphs

/home/tianye/code/Script/GMT/TXT2CPT_vel_4decimal $event'_am_laplace.txt.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_sta $cpt_path$event'_am_laplace.txt.HD.cpt' $event'_am_laplace.txt.HD' $event'_am.txt_v2' $argv[4]
mv $event'_am_laplace.txt.HD.ps' $event.graphs
exit
rm -f 'vel_crctd_'$event'.txt.HD'
foreach am_file (`awk '{print}' $event'_am_laplace.txt.HD' | sed s/' '/'@'/g`)
set long=`echo $am_file | cut -d@ -f1`
set lati=`echo $am_file | cut -d@ -f2`
if(! `awk -v long=$long -v lati=$lati '(long-$1)^2+(lati-$2)^2<0.00001' 'slow_term_'$event'.txt.HD' | wc -l`)then
continue
endif
set am_term=`echo $am_file | cut -d@ -f3`
echo $am_term
set slow_term=`awk -v long=$long -v lati=$lati '(long-$1)^2+(lati-$2)^2<0.00001 {print $3}' 'slow_term_'$event'.txt.HD'`
#echo $am_term $slow_term
set vel_crctd=`echo $slow_term $am_term | awk '{printf "%.3f\n",sqrt(1/($1-$2))}'`
echo $long $lati $vel_crctd >> 'vel_crctd_'$event'.txt.HD'
end
/home/tianye/code/Script/GMT/TXT2CPT_vel 'vel_crctd_'$event'.txt.HD'
/home/tianye/code/Script/GMT/C_plot_kernel_input_sta $cpt_path'vel_crctd_'$event'.txt.HD.cpt' 'vel_crctd_'$event'.txt.HD' $event'_am.txt_v2' $argv[4]
mv 'vel_crctd_'$event'.txt.HD.ps' $event.graphs
