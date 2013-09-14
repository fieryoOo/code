#!/bin/csh
if ($#argv != 2)then
echo "usage: Pick_data [sac_path] [v1 or v2]"
exit
endif

cd $argv[1]
foreach period (22 30 40 50 60 70 80 90 100)
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis/amp_diff_'$period's_'$argv[2]
echo "working on "$period"sec period..."
rm -f sta_avg_std_table
rm -f sta_avg
rm -f sta_avg2
rm -f sta_std
rm -f sta_std2
foreach file (`ls ????_pctg`)
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
echo $avg >> sta_avg
echo $avg2 >> sta_avg2
echo $std >> sta_std
echo $std2 >> sta_std2
#exit
end
#exit
cd ../..
end
