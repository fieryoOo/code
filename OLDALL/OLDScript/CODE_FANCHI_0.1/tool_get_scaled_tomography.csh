#!/bin/csh
#echo on
if ($#argv != 1) then
  echo "USAGE: tool_do_AtoD  period"
  exit 1
endif

set per=$argv[1]

@ dis = 12 * $per
cd $per"sec_10snr_"$dis"dis"
#foreach scale ( 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)
foreach scale ( 1.0 )
#/home/linf/SCRIPTS/bin_64/slow_laplace_maps_to_iso_map_ani_data_v5_n50_270_noamp_scale event_lst_cor_30_100 0 360 18 $per"s_iso_ani_v50" $scale
#/home/weisen/PROGS_64/EIKONAL_earthquake/slow_laplace_maps_to_iso_map_ani_data_v5_n50_270_noamp_scale event_lst_cor_30_100 0 360 18 $per"s_iso_ani_v50" $scale
/home/weisen/PROGS_64/EIKONAL_earthquake/GRID_1/slow_laplace_maps_to_iso_map_ani_data_v5_n50_270_noamp_scale_input_small_region_cv1 event_lst_cor_30_100 0 360 18 $per"s_iso_ani_v50" $scale 248 34 102 92
end
cd ..

