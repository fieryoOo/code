#!/bin/csh
#echo on
if ($#argv != 1) then
  echo "USAGE: tool_do_AtoD  period"
  exit 1
endif

set per = $1
#(  24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 58 60 64 68 72 76 80 84 88 92 96 100)
@ dis = 12 * $per
cd $per"sec_10snr_"$dis"dis"

ls *ph.txt -lh | awk '{if ($5>0) print $8}' | awk -F. '{print $1}' > good_event.dat

foreach event ( `awk '{print $1}' good_event.dat` )


#/home/linf/SCRIPTS/bin_64/lf_correct_2pi_v1 $event".ph.txt" $per
#/home/weisen/PROGS_64/EIKONAL_earthquake/lf_correct_2pi_v1 $event".ph.txt" $per
/home/weisen/PROGS_64/EIKONAL_earthquake/GRID_1/lf_correct_2pi_v1_small_region_cv1 $event".ph.txt" $per
../C_plot_travel $event".ph.txt_v1"
../C_plot_travel_am  $event".ph.txt_v1"
#/home/linf/SCRIPTS/bin_64/correct_travel_time_curvature_v1_270 $event $per
#/home/weisen/PROGS_64/EIKONAL_earthquake/correct_travel_time_curvature_v1_270 $event $per 
/home/weisen/PROGS_64/EIKONAL_earthquake/GRID_1/correct_travel_time_curvature_v1_270_input_small_region_cv2 $event $per 248 34 102 92   # 34 43 
#../C_plot_travel $event".ph.txt_v2"
#../C_plot_travel_T0.2  $event".ph.txt_v2"
#../C_plot_travel_am  $event".ph.txt_v2"
end
#ls 20*.ph.txt_v2 | cut -d. -f1 > event_lst
#/home/linf/SCRIPTS/bin_64/travel_time_to_slow_map_v4_v3_v2to10_270 event_lst $per
#/home/linf/SCRIPTS/bin_64/slow_maps_to_iso_map_ani_data_v5_n50_270 event_lst 0 360 18 $per"s_iso_ani_v1"
#rm 20*HD*
cd ..
