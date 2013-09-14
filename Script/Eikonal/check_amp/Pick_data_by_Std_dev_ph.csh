#!/bin/csh
if ($#argv != 1)then
echo "usage: Pick_data [sac_path]"
exit
endif

cd $argv[1]
foreach period (22 30 40 50 60 70 80 90 100)
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis/ph_diff_'$period's/'
echo "working on "$period"sec period..."
rm -f sta_avg_std_table
foreach file (`ls ????`)
#set file='SC69_pctg'
set avg=`awk 'BEGIN{a=0}{a=a+$1}END{print a/NR}' $file`
set std=`awk -v a=$avg 'BEGIN{s=0}{s=s+($1-a)^2}END{print sqrt(s/(NR-1))}' $file`
set b=`echo "$avg-$std*1.8" | bc`
set e=`echo "$avg+$std*1.8" | bc`
#awk -v b=$b -v e=$e '{if($1<b || $1>e){print $1}}' $file | wc -l
awk -v b=$b -v e=$e '{if($1>b && $1<e){print $1}}' $file > $file'.picked'
if (`more $file'.picked' | wc -l` < 10)then
continue
endif
set avg2=`awk 'BEGIN{a=0}{a=a+$1}END{print a/NR}' $file'.picked'`
set std2=`awk -v a=$avg2 'BEGIN{s=0}{s=s+($1-a)^2}END{print sqrt(s/(NR-1))}' $file'.picked'`
echo $file' '$avg' '$std' '$avg2' '$std2 >> sta_avg_std_table
#exit
end
#exit
cd ../..
end
