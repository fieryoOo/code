#!/bin/csh
if ( $#argv != 1)then
echo "USAGE: get_event [SAC_path]"
exit
endif
cd $argv[1]
foreach per ( 22 30 40 50 60 70 80 90 100 )
@ dis = 12 * $per
/home/tianye/code/Programs/Earthquake/get_event_period_trvt_amp $per sta_dist.lst 10 $dis $argv[1] > test1.log
end
