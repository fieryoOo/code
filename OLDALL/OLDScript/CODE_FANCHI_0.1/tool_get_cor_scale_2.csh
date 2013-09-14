#!/bin/csh
#echo on
if ($#argv != 1) then
  echo "USAGE: tool_do_AtoD  period"
  exit 1
endif

set per=$argv[1]

@ dis = 12 * $per
cd $per"sec_10snr_"$dis"dis"
awk '{print $3}' log.log > log_cor.log
awk '{print 1/$4}' log.log > log_scale.log
#/home/linf/SCRIPTS/bin_64/get_hist_percent_v1 log_cor.log -1 1 40 > log_cor.hist
#/home/linf/SCRIPTS/bin_64/get_hist_percent_v1 log_scale.log 0 80 160 > log_scale.hist
#/home/linf/SCRIPTS/bin_64/get_hist_percent_v1_midpoint log_cor.log -1 1 40 > log_cor.mid
#/home/linf/SCRIPTS/bin_64/get_hist_percent_v1_midpoint log_scale.log 0 80 160 > log_scale.mid
/home/weisen/PROGS_64/EIKONAL_earthquake/get_hist_percent_v1 log_cor.log -1 1 40 > log_cor.hist
/home/weisen/PROGS_64/EIKONAL_earthquake/get_hist_percent_v1 log_scale.log 0 80 160 > log_scale.hist
/home/weisen/PROGS_64/EIKONAL_earthquake/get_hist_percent_v1_midpoint log_cor.log -1 1 40 > log_cor.mid
/home/weisen/PROGS_64/EIKONAL_earthquake/get_hist_percent_v1_midpoint log_scale.log 0 80 160 > log_scale.mid
awk '{if($3>0.3) print $1}' log.log | cut -d_ -f2 > event_lst_cor_30_100
cd ..

