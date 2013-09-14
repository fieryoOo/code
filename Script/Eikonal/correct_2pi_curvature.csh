#!/bin/csh
#echo on
if ($#argv != 2) then
  echo "USAGE: tool_do_AtoD [region_infile] [sac_path]"
  exit 1
endif

cd $argv[2]

set REG = `more $argv[1]`
#echo $REG

set flag=0
#foreach per ( 10 15 16 20 25 )
foreach per ( 10 12 15 20 25 30 40 50 70 90 )
@ dis = 4 * $per
cd $per"sec_8snr_"$dis"dis"
if( $flag == 0 )then
foreach event(`ls 2007*.ph.txt | cut -d. -f1`)
if( ! -e /home/tianye/Model/PHASE_map/Event_Pre/$event )then
/home/tianye/code/Programs/Model/Event_Modle.csh ../station.lst ../event.loc
endif
end
set flag=1
endif
foreach event(`ls 2007*.ph.txt | cut -d. -f1`)
/home/tianye/code/Programs/Earthquake/correct_2pi_model $event".ph.txt" $per
if( $per == 25)then
grep -v 'B18A' $event".ph.txt_v1" > temp1
grep -v 'I15A' temp1 > $event".ph.txt_v1"
rm -f temp1
endif
set lonl=`more $argv[1] | cut -d/ -f1 | sed s/'-R'/''/`
set lonh=`more $argv[1] | cut -d/ -f2`
set latl=`more $argv[1] | cut -d/ -f3`
set lath=`more $argv[1] | cut -d/ -f4`
awk -v lonl=$lonl -v lonh=$lonh -v latl=$latl -v lath=$lath '$1>lonl && $1<lonh && $2>latl && $2<lath' $event".ph.txt_v1" > temp
mv temp $event".ph.txt_v1"
/home/tianye/code/Script/GMT/C_plot_travel_positive $event".ph.txt_v1" $argv[1] 0.1
/home/tianye/code/Script/GMT/C_plot_travel_am  $event".ph.txt_v1" $argv[1]
#/home/tianye/code/Programs/Earthquake/correct_travel_time_curvature_v1_270 $event $per 50 0.02 0.05

end
exit

rm -f `wc $event'.ph.txt_v2' | grep "0       0       0" | awk '{print $4}'`
if( ! -e $event'.ph.txt_v2' ) continue
/home/tianye/code/Script/GMT/C_plot_travel_positive $event'.ph.txt_v2' $argv[1] 0.1
/home/tianye/code/Script/GMT/C_plot_travel_T0.2 $event'.ph.txt_v2' $argv[1] 0.1
/home/tianye/code/Script/GMT/C_plot_travel_am $event'.ph.txt_v2' $argv[1]
end
cd ..
end
