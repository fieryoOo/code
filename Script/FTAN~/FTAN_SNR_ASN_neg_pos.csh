#!/bin/csh
if ( $#argv != 2)then
echo "USAGE: FTAN [SAC_path] [COR or DCV]"
exit
endif

set tp = $argv[2]
cd $argv[1]
#foreach dir (`ls -d */`)
foreach dir (`awk '{print $1}' station.lst`)
#foreach dir (T25A)
if ( ! -e $dir )continue
cd $dir
if ( `ls $tp*'.SAC'| wc -l` == 0 )then
cd ..
continue
endif
rm -f $tp*'.SAC_'*'_DISP.'?
rm -f $tp*'.SAC_'*'amp_snr'
rm -f $tp*'.SAC_'*'_AMP'
ls $tp*'.SAC' > sac.lst
#touch_sac sac.lst
set param = '-1 0.6 4.5 1 40 20 1. 1.0 0.2 1.'
awk -v param="$param" '{print param,$1,"1" }' sac.lst > param_R.dat
rm -f sac.lst
/home/tianye/code/Programs/FTAN_amp_snr_1Dmodel/aftani_c_pgl_amp param_R.dat /home/tianye/code/Programs/head/OBS_phvel.dat
#/home/tianye/code/Programs/FTAN_amp_snr_1Dmodel/aftani_c_pgl_amp param_R.dat /home/tianye/code/Programs/head/overtone.dat
#foreach file ( `ls *.SAC_2_DISP.1` )
#if( `awk '$5>5.0' $file | wc -l` > 0 )then
#set file2=`echo $file | sed s/'_2_DISP.0'/''/g`
#rm $file2'_'*
#endif
#end
foreach file ( `ls $tp'_'*'_amp_snr'` )
if( `awk '$1<23 && $1>6' $file | awk -v file=$file 'BEGIN{flag=1}{if($3>5 || $5>5){flag=0}}END{print flag}'` )then
set file2=`echo $file | sed s/'_amp_snr'/''/g`
rm -f $file2'_'*
endif
end
### redo FTAN on both sides
foreach file ( `ls *.SAC_2_DISP.1` )
#awk '{print $3,$5}' $file > temp.dat
set sac=`echo $file | sed s/'SAC_2_DISP.1'/'SAC'/g`
set sacn=$sac'_neg'
set sacp=$sac'_pos'
sac << END
cut -3000 0
r $sac
reverse
ch b 0
w $sacn
cut 0 3000
r $sac
w $sacp
quit
END
ls $sac'_neg' > sac.lst
ls $sac'_pos' >> sac.lst
#touch_sac sac.lst
awk -v param="$param" '{print param,$1,"0" }' sac.lst > param_R.dat
/home/tianye/code/Programs/FTAN_amp_snr_1Dmodel/aftani_c_pgl_amp param_R.dat /home/tianye/code/Programs/head/OBS_phvel.dat
end
cd ..
end
