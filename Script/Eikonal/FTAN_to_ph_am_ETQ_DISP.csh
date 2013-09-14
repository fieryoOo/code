#!/bin/csh
if ( $#argv != 2)then
echo "USAGE: get_event [SAC_path] [sta_dist.lst]"
exit
endif
cd $argv[1]
foreach per (60)
#foreach per ( 8 10 14 20 30 )

@ dis = 4 * $per
/home/tianye/code/Programs/Earthquake/get_event_period_trvt_amp_DISP $per $argv[2] 6 $dis $argv[1]

foreach file (`ls $per'sec_6snr_'$dis'dis/20'*'.ph.txt'`)
if( `more $file | wc -l` < 80 )then
rm -f $file
continue
endif
grep -v 'NaN' $file > temp
set ev=`echo $file | cut -d/ -f2 | cut -d. -f1`
if(! `grep -c $ev event.loc`)then
rm -f $file
continue
endif
set loc=`grep -m1 $ev event.loc | awk '{print $2,$3}'`
echo $loc 0.0 NaN NaN $ev > $file
cat temp >> $file
end

end
rm -f temp
