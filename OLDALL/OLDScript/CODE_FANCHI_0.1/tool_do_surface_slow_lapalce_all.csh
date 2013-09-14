#!/bin/csh
#echo on
if ($#argv != 1) then
  echo "USAGE: tool_do_AtoD  period"
  exit 1
endif

set per=$argv[1]

@ dis = 12 * $per
cd $per"sec_10snr_"$dis"dis"
 wc 20*ph.txt_v2 | grep "0       0       0" | awk '{print "rm",$4}' > rm_0_0_0.csh
csh rm_0_0_0.csh
ls 20*ph.txt_v2 | awk '{print "../C_plot_travel",$1,"\n../C_plot_travel_T0.1",$1}' > plot_all.csh
csh plot_all.csh
ls 20*.ph.txt_v2 | cut -d. -f1 > event_lst
#/home/linf/SCRIPTS/bin_64/travel_time_to_slow_map_v4_v3_v2to10_270_noampHD event_lst $per
#/home/linf/SCRIPTS/bin_64/slow_maps_to_iso_map_ani_data_v5_n50_270_noamp event_lst 0 360 18 $per"s_iso_ani_v1"
/home/weisen/PROGS_64/EIKONAL_earthquake/GRID_1/travel_time_to_slow_map_v4_v3_v2to10_270_noampHD_input_small_region_cv1 event_lst $per 248 34 102 92
/home/weisen/PROGS_64/EIKONAL_earthquake/GRID_1/slow_maps_to_iso_map_ani_data_v5_n50_270_noamp_input_small_region_cv1 event_lst 0 360 18 $per"s_iso_ani_v1" 248 34 102 92


grep -v " 0 999" $per"s_iso_ani_v1.iso"  | awk '{print $1,$2,$3}' > $per"s_iso_ani_v1.1"
cd ..

