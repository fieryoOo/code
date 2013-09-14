#!/bin/csh
if ( $#argv != 1 )then
echo "usage: Calc_overall.csh [sac_path]"
exit
endif
cd $argv[1]
rm -f ph_diff_avg_period_table
foreach period (22 30 40 50 60 70 80 90 100)
@ dis = $period * 12
set avg=`awk 'BEGIN{a=0}{a=a+$4}END{print a/NR}' $period'sec_10snr_'$dis'dis/ph_diff_'$period's/sta_avg_std_table'`
set std=`awk 'BEGIN{a=0}{a=a+$5}END{print a/NR}' $period'sec_10snr_'$dis'dis/ph_diff_'$period's/sta_avg_std_table'`
echo $period's '$avg' '$std >> ph_diff_avg_period_table
end
