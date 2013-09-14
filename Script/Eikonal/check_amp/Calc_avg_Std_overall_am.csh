#!/bin/csh
if ( $#argv != 2 )then
echo "usage: Calc_overall.csh [sac_path] [v1 or v2]"
exit
endif
cd $argv[1]
rm -f amp_diff_avg_period_table_$argv[2]
foreach period (22 30 40 50 60 70 80 90 100)
@ dis = $period * 12
set avg=`awk 'BEGIN{a=0}{a=a+$4}END{print a/NR}' $period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[2]'/sta_avg_std_table'`
set std=`awk 'BEGIN{a=0}{a=a+$5}END{print a/NR}' $period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[2]'/sta_avg_std_table'`
echo $period's '$avg' '$std >> amp_diff_avg_period_table_$argv[2]
end
