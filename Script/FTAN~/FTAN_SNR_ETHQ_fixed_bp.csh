#!/bin/csh
if ( $#argv != 1)then
echo "USAGE: FTAN [SAC_path]"
exit
endif

cd $argv[1]
#foreach dir (`ls -d */`)

rm -f moevent.loc
foreach evloc (`awk '{print $1"@"$2"@"$3}' event.loc`)
set ev=`echo $evloc | cut -d@ -f1`
set lon=`echo $evloc | cut -d@ -f2`
set lat=`echo $evloc | cut -d@ -f3`
if ( ! -e $ev )continue
cd $ev
echo "Removing old disp and snr files for event "$ev"..."
rm -f $ev*'.sac_2_DISP.'?
rm -f $ev*'.sac_amp_snr'
if ( `ls $ev*'.sac'| wc -l` == 0 )then
cd ..
continue
endif
cd ..
if( ! -e /home/tianye/Model/PHASE_map/Event_Pre/$ev )then
echo $ev $lon $lat >> moevent.loc
endif
end

echo "Computing phc modle..."
if( -e moevent.loc )/home/tianye/code/Programs/Model/Event_Modle.csh ~/station_TA.lst moevent.loc
#rm -f moevent.loc

foreach ev (`awk '{print $1}' event.loc`)
cd $ev
echo "Working on event "$ev"..."
ls $ev*'.sac' > sac.lst
#touch_sac sac.lst
awk '{print "0 2.0 5.0 5 100 20 1 0.5 0.2 2 "$1" 0" }' sac.lst > param_R.dat
rm -f sac.lst

#ls 'COR'*'.SAC' | awk '{print "-1 2.0 5.0 5 50 20 1 0.5 0.2 2",$1,"1" }' > param_R.dat
/home/tianye/code/Programs/FTAN_amp_snr_ethq_fixed_bp/aftani_c_pgl_amp param_R.dat /home/tianye/Model/PHASE_map/Event_Pre/
foreach file ( `ls $ev*'.sac_2_DISP.0'` )
if( `awk '$5>4.5' $file | wc -l` > 0 )then
set file2=`echo $file | sed s/'_2_DISP.0'/''/g`
rm $file2'_'*
endif
end
foreach file ( `ls $ev*'_amp_snr'` )
if( `awk '$1<23 && $1>12' $file | awk -v file=$file 'BEGIN{flag=1}{if($3>5 || $5>5){flag=0}}END{print flag}'` )then
set file2=`echo $file | sed s/'_amp_snr'/''/g`
rm $file2'_'*
endif
end
cd ..
end
