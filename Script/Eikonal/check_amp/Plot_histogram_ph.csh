#!/bin/csh
if ($#argv != 1)then
echo "usage: Plot_hist [sac_path]"
exit
endif

cd $argv[1]
set net=`echo $argv[1] | awk -F/ '{print $(NF-1)}'`
foreach period (22 30 40 50 60 70 80 90 100)
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis/ph_diff_'$period's/'
cp /home/tianye/code/Script/GMT/gmtdefaults4 .gmtdefaults4

set avg=`awk 'BEGIN{a=0}{a+=$4}END{print a/NR}' sta_avg_std_table`
awk '{print $4}' sta_avg_std_table | pshistogram -B0.2/5:."$net  mean $avg": -R-.6/.6/0/30 -Z1 -Jx19.5c/0.5c -W.1 -L1.0p,red -F -Y2.5 -X4.1 > 'sta_avg_std_table.ps'

foreach file (`ls ????.picked`)
set sta=`echo $file | cut -d. -f1`
set temp=`awk -v sta=$sta '$1==sta {print $4"@"$5}' sta_avg_std_table`
set avg=`echo $temp | cut -d@ -f1`
set std=`echo $temp | cut -d@ -f2`
pshistogram $file -B1/5:."$period sec  $file  mean $avg  std $std": -R-3/3/0/30 -Z1 -Jx3.9c/0.5c -W0.2 -L1.0p,red -F -Y2.5 -X4.1 > $file'.ps'
end
cd ../..
end
