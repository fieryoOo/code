#!/bin/csh
if ( $#argv != 1)then
echo "USAGE: FTAN [SAC_path]"
exit
endif

cd $argv[1]
foreach dir (`ls -d 20*/`)
cd $dir
ls '20'*'.sac' > sac.lst
touch_sac sac.lst
ls '20'*'.sac' | awk '{print "0 2.0 5.0 5 100 20 1 0.5 0.2 2 "$1" 0" }' > param_R.dat
/home/tianye/code/Programs/FTAN_amp_snr_nomod/aftani_c_pgl_amp param_R.dat
#ls '20'*'.sac' | awk '{print "0 2.0 5.0 5 150 20 1 0.5 0.2 2 "$1 }' > param_R.dat
#/home/tianye/code/Programs/FTAN_amp/aftani_c_pgl_amp param_R.dat
#/home/tianye/code/Programs/SPEC_SNR_RMS/lf_spec_snr_rms_fast sac.lst
rm -f sac.lst

cd ..
end
