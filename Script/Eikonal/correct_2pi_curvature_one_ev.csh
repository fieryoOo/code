#!/bin/csh
#echo on
if ($#argv != 2) then
  echo "USAGE: tool_do_AtoD [region_infile] [sac_path]"
  exit 1
endif

cd $argv[2]

set REG = `more $argv[1]`
#echo $REG
foreach per ( 80 )
@ dis = 12 * $per
cd $per"sec_10snr_"$dis"dis"

#foreach event(`ls 20*.ph.txt | cut -d. -f1`)
set event = 20091027000446
#/home/tianye/code/Programs/Earthquake/correct_2pi_v1 $event".ph.txt" $per 20
#/home/tianye/code/Script/GMT/C_plot_travel $event".ph.txt_v1" $argv[1]
#/home/tianye/code/Script/GMT/C_plot_travel_am  $event".ph.txt_v1" $argv[1]
/home/tianye/code/Programs/Earthquake/correct_travel_time_curvature_v1_270 $event $per 20 0.02 0.04
#end

rm -f `wc $event'.ph.txt_v2' | grep "0       0       0" | awk '{print $4}'`
/home/tianye/code/Script/GMT/C_plot_travel $event'.ph.txt_v2' $argv[1]
/home/tianye/code/Script/GMT/C_plot_travel_T0.2 $event'.ph.txt_v2' $argv[1]
/home/tianye/code/Script/GMT/C_plot_travel_am $event'.ph.txt_v2' $argv[1]

cd ..
end
