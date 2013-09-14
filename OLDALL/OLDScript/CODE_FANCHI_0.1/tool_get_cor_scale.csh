#!/bin/csh
#echo on
if ($#argv != 1) then
  echo "USAGE: tool_do_AtoD  period"
  exit 1
endif

set per=$argv[1]

@ dis = 12 * $per
cd $per"sec_10snr_"$dis"dis"
#/home/linf/SCRIPTS/src_32_static/Cal_Laplace/fit_tr_amp_to_second_order_v50_fit_scale/fit_tr_amp_to_second_order event_lst $per"s_iso_ani_v1.1" > log.log
#/home/weisen/PROGS_64/EIKONAL_earthquake/fit_tr_amp_to_second_order_v50_fit_scale/fit_tr_amp_to_second_order event_lst $per"s_iso_ani_v1.1" > log.log
/home/weisen/PROGS_64/EIKONAL_earthquake/GRID_1/fit_tr_amp_to_second_order_v50_fit_scale/fit_tr_amp_to_second_order_input_small_region event_lst $per"s_iso_ani_v1.1" 248 34 102 92 > log.log
cd ..

