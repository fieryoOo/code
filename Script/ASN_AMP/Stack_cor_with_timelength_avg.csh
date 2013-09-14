#!/bin/csh
#if(#$argv != 2)then
#echo "Usage: Calc_cor_time [sta_1] [sta_2]"
#exit
#endif
#cd /utera/tianye/data_avg_monthly_stack
#foreach info (`awk '{print $1"@"$2"@"}' H19A_pair.lst`)
#set sta1=`echo $info | cut -d@ -f1`
#set sta2=`echo $info | cut -d@ -f2`
#set dir=`echo $sta1"_"$sta2`
#echo "Stacking station pair "$dir
#cd $dir
#set flag=0
#foreach year (2008 2009 2010)
#foreach month (JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC)
#set file=$year'_'$month'_COR_'$dir'.SAC'
#if(! -e $file)continue
#set day_num=`grep $month rec/$dir'_'$year'.txt' | awk '{print $4}'`
##echo $year $month $day_num
#if($flag == 0)then
##echo $day_num
#cp $file $dir'_stacked_'$day_num'_days.SAC'
#set day_sum=$day_num
#set flag=1
#else 
#set stack_name=`ls $dir'_stacked_'$day_sum'_days.SAC'`
#set day_sum=`echo $day_num $day_sum | awk '{print $1+$2}'`
#set new_name=$dir'_stacked_'$day_sum'_days.SAC'
##echo $file
##echo $stack_name
#sac << END
#BINOPERR NPTS WARNING
#cut -3000 3000
#r $stack_name
#addf $file
#cut off
#w $new_name
#quit
#END
#endif
#end
#end
#cd ..
#end

cd /utera/tianye/data_avg_monthly_stack
rm -f stacked_sac.lst
ls -d */ | grep 'LRL' > stack_pair.lst
ls -d */ | grep 'H19A' >> stack_pair.lst
ls -d */ | grep 'X25A' >> stack_pair.lst
foreach dir (`more stack_pair.lst`)
ls $dir*'_stacked_'*'days.SAC' >> stacked_sac.lst
end
#/home/tianye/code/Programs/seed2sac_earthquake/lf_touch_sac stacked_sac.lst

/home/tianye/code/Programs/CORRECT_ASN_AMP/nosnr_nonoise/amp_length_normalize stacked_sac.lst

foreach info (`awk '{print $1"@"$2}' LRL_H19A_X25A_pair.lst`)
set sta1=`echo $info | cut -d@ -f1`
set sta2=`echo $info | cut -d@ -f2`
set dir=`echo $sta1"_"$sta2`
#set dist=`echo $info | awk '{print $3}'`
echo "Extracting amp_noise_days info..Working on station pair "$dir'...'
cd $dir
set file_1=`ls $dir'_stacked_'*'_days.SAC_amp_l_norm.txt' | head -1`
foreach period (`awk '{print $1}' $file_1`)
rm -f $period'_noisesqr_days_noise_ampm_ampc.txt'
foreach stack (`ls $dir'_stacked_'*'_days.SAC_amp_l_norm.txt'`)
set time_l=`echo $stack | cut -d_ -f4`
set data=`awk -v per=$period '$1==per' $stack`
set amp_max=`echo $data | awk '{print $3}'`
set amp_crtd=`echo $data | awk '{print $2}'`
set noise_sqr=`echo $data | awk '{print $4**2}'`
set noise_rms=`echo $data | awk '{print $4}'`
echo $noise_sqr $time_l $noise_rms $amp_max $amp_crtd >> $period'_noisesqr_days_noise_ampm_ampc.txt'
end
end
cd ..
end

rm -f combined.lst
foreach sta (LRL H19A X25A)
echo "Combining data for station "$sta"..."
set dir1=`ls -d */ | grep -m1 $sta | cut -d/ -f1`
set file1=`ls $dir1'/'$dir1'_stacked_'*'_days.SAC_amp_l_norm.txt' | head -1`
rm -f periods.txt
foreach period (`awk '{print $1}' $file1`)
echo $period >> periods.txt
rm -f combined_$sta'_'$period'_noisesqr_days_noise_ampm_ampc.txt'
foreach pair (`ls -d */ | grep $sta`)
more $pair'/'$period'_noisesqr_days_noise_ampm_ampc.txt' >> combined_$sta'_'$period'_noisesqr_days_noise_ampm_ampc.txt'
end
echo combined_$sta'_'$period'_noisesqr_days_noise_ampm_ampc.txt' >> combined.lst
end
end

/home/tianye/code/Programs/CORRECT_ASN_AMP/calc_avg_daynum_range 60 combined.lst

foreach file ( `ls -d *combined_*_noisesqr_days_noise_ampm_ampc.txt` )
set file_new=`echo $file | sed s/'ampc'/'ampc_2nd_year'/g`
awk '$2>80 && $2<650' $file > $file_new
end
