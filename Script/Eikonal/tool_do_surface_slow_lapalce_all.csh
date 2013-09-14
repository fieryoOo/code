#!/bin/csh
if ( $#argv != 2 )then
echo "usage: do_laplace [region_inputfile] [sac_path]"
exit
endif

cd $argv[2]

foreach per ( 22 30 40 50 60 70 80 90 100 )

echo "\n\nworking on "$per"sec period...\n"
@ dis = 12 * $per
cd $per"sec_10snr_"$dis"dis"
wc 20*ph.txt_v2 | grep "0       0       0" | awk '{print "rm",$4}' > rm_0_0_0.csh
echo "removing 0_0_0 files..."
csh rm_0_0_0.csh
ls 20*ph.txt_v2 | awk -v region=$argv[1] '{print "/home/tianye/code/Script/GMT/C_plot_travel",$1,region,"\n/home/tianye/code/Script/GMT/C_plot_travel_T0.2",$1,region}' > plot_all.csh
echo "converting ph to HD..."
csh plot_all.csh
ls 20*.ph.txt_v2 | cut -d. -f1 > event_lst
#/home/tianye/code/Programs/Earthquake/travel_time_to_slow_map_v4_v3_v2to10_270_noampHD event_lst $per
#/home/tianye/code/Programs/Earthquake/slow_maps_to_iso_map_ani_data_v5_n50_270_noamp event_lst 0 360 18 $per"s_iso_ani_v1"
#grep -v " 0 999" $per"s_iso_ani_v1.iso"  | awk '{print $1,$2,$3}' > $per"s_iso_ani_v1.1"

ls 20*ph.txt_v2 | awk -v region=$argv[1] '{print "/home/tianye/code/Script/GMT/C_plot_travel_am",$1,region}' > plot_all_am.csh
echo "converting am to HD..."
csh plot_all_am.csh
foreach event (`more event_lst`)
rm -f $event'_am_laplace.txt.HD'
end
/home/tianye/code/Programs/Earthquake/travel_time_to_slow_map_v4_v3_v2to10_270_noampHD event_lst $per
/home/tianye/code/Programs/Earthquake/slow_maps_to_iso_map_ani_data_v5_n50_270_noamp event_lst 0 360 18 $per"s_iso_ani_v1"
/home/tianye/code/Programs/Earthquake/amp_HD_to_amp_gradient_HD_to_amp_laplace_HD event_lst $per
/home/tianye/code/Programs/Earthquake/slow_laplace_maps_to_iso_map_ani_data_v5_n50_270_noamp_scale event_lst 0 360 18 $per"s_laplace_iso_ani_v1" 1
grep -v '0 9999 0' $per's_iso_ani_v1.iso' | awk '{print $1,$2,$3}' > $per's_iso_ani_v1.1'
grep -v " 0 9999" $per"s_laplace_iso_ani_v1_scale_1.iso"  | awk '{print $1,$2,$3}' > $per"s_laplace_iso_ani_v1.1_scale_1.iso"
cd ..
#exit
end
