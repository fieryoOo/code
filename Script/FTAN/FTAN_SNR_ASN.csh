#!/bin/csh
if ( $#argv != 1)then
echo "USAGE: FTAN [SAC_path]"
exit
endif

cd $argv[1]
foreach dir ( P22A )
#foreach dir (`awk '{print $1}' station.lst`)
if ( ! -e $dir )continue
cd $dir
rm -f 'COR'*'.SAC_'?'_DISP.'?
rm -f 'COR'*'.SAC_amp_snr'
rm -f 'COR'*'.SAC_'?'_AMP'
if ( `ls 'COR'*'.SAC'| wc -l` == 0 )then
cd ..
continue
endif
ls 'COR'*'.SAC' > sac.lst
touch_sac sac.lst
rm -f sac.lst
ls 'COR_'*'.SAC' | awk '{print "-1 2.0 5.0 5 40 20 1. 0.5 0.2 2",$1,"1" }' > param_R.dat
/home/tianye/code/Programs/FTAN_amp_snr/aftani_c_pgl_amp param_R.dat /home/tianye/Model/PHASE_map/US_intersta_Pre/
###
cd ..
continue
###
foreach file ( `ls *.SAC_2_DISP.1` )
if( `awk '$5>5.0' $file | wc -l` > 0 )then
set file2=`echo $file | sed s/'_2_DISP.0'/''/g`
rm $file2'_'*
endif
end
foreach file ( `ls COR_*_amp_snr` )
if( `awk '$1<23 && $1>6' $file | awk -v file=$file 'BEGIN{flag=1}{if($3>5 || $5>5){flag=0}}END{print flag}'` )then
set file2=`echo $file | sed s/'_amp_snr'/''/g`
rm $file2'_'*
endif
end
cd ..
end
