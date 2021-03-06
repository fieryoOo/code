#!/bin/csh
if ( $#argv != 3 )then
echo "usage: Avg_l_normd_Amp.csh [sac_path] [sta_p_dist_days_table] [sta_lst]"
exit
endif

cd $argv[1]
set label=`echo $argv[3] | awk -F/ '{print $NF}' | cut -d. -f1 | cut -d_ -f2`
#rm -f sta_for_amp_norm_days.lst
#rm -f sta_for_amp_norm.lst
#set N=`more $argv[2] | wc -l`
#set i=0
#set temp2=0
#foreach sta_p (`awk '{print $1"@"$2"@"$4}' $argv[2]`)
#set temp=`echo $i $N | awk '{printf "%.0f",$1*100/$2}'`
#if ($temp != $temp2 )then
#echo "Producing norm file list..."$temp2" percent finished..."
#endif
#set temp2=$temp
#
#set sta1=`echo $sta_p | cut -d@ -f1`
#set sta2=`echo $sta_p | cut -d@ -f2`
#set dayn=`echo $sta_p | cut -d@ -f3`
#if ( -e $sta1'/COR_'$sta1'_'$sta2.SAC )then
#echo $sta1'/COR_'$sta1'_'$sta2.SAC $dayn >> sta_for_amp_norm_days.lst
#else if ( -e $sta2'/COR_'$sta2'_'$sta1.SAC )then
#echo $sta2'/COR_'$sta2'_'$sta1.SAC $dayn >> sta_for_amp_norm_days.lst
#endif
#@ i += 1
#end
#
#echo "touching those fucking foolish SAC files..."
#awk '{print $1}' sta_for_amp_norm_days.lst > sta_for_amp_norm.lst
#/home/tianye/code/Programs/seed2sac_earthquake/lf_touch_sac sta_for_amp_norm.lst
echo "Normalizing the amplitude for each station pairs..."
/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/amp_length_normalize sta_for_amp_norm_days.lst

mkdir -p $label'_normd_amp'
foreach sta_info (`awk '{print $1"@"$2"@"$3}' $argv[3]`)
set sta=`echo $sta_info | cut -d@ -f1`
set long=`echo $sta_info | cut -d@ -f2`
set lati=`echo $sta_info | cut -d@ -f3`
rm -f $label'_normd_amp/'$sta'_avgd_daynormd_amp'
rm -f $label'_normd_amp/'$sta'_avgd_noisenormd_amp'
rm -f $sta'_file.lst'
foreach file ( `grep '_'$sta sta_for_amp_norm.lst | awk '{print $1"_amp_l_norm.txt"}'` )
if(-e $file)echo $file >> $sta'_file.lst'
end
if ( `more $sta'_file.lst' | wc -l` < 5 )then
echo $sta": less than 5 stations!"
rm -f $sta'_file.lst'
continue
endif
echo "Calculating the average amplitude...Working on station "$sta"..."
set file1=`head -1 $sta'_file.lst'`
set dis_sta1=`echo $file1 | cut -d_ -f2`
set dis_sta2=`echo $file1 | cut -d_ -f3 | cut -d. -f1`
set dist=`grep $dis_sta1 $argv[2] | grep -m1 $dis_sta2 | awk '{print $3}'`
set weight_sum=`echo $dist | awk '{print 1/$1**2}'`
awk -v dist=$dist '{print $1,$2/dist**2}' $file1 > add_temp_day
awk -v dist=$dist '{print $1,$3/dist**2}' $file1 > add_temp_noise
set i=0
foreach file (`more $sta'_file.lst'`)
if(i == 0)then
@ i += 1
continue
endif
set dis_sta1=`echo $file | cut -d_ -f2`
set dis_sta2=`echo $file | cut -d_ -f3 | cut -d. -f1`
set dist=`grep $dis_sta1 $argv[2] | grep -m1 $dis_sta2 | awk '{print $3}'`
set weight_sum=`echo $dist $weight_sum | awk '{print 1/$1**2+$2}'`
awk -v file=$file -v dist=$dist '{a=$2;getline < file;print $1,$2/dist**2+a}' add_temp_day > add_temp_day_2
awk -v file=$file -v dist=$dist '{a=$2;getline < file;print $1,$3/dist**2+a}' add_temp_noise > add_temp_noise_2
mv add_temp_day_2 add_temp_day
mv add_temp_noise_2 add_temp_noise
@ i += 1
end
rm -f $sta'_file.lst'
echo $i
echo $sta $long $lati > $label'_normd_amp/'$sta'_avgd_daynormd_amp'
echo $sta $long $lati > $label'_normd_amp/'$sta'_avgd_noisenormd_amp'
awk -v weight=$weight_sum '{print $1,$2/weight}' add_temp_day >> $label'_normd_amp/'$sta'_avgd_daynormd_amp'
awk -v weight=$weight_sum '{print $1,$2/weight}' add_temp_noise >> $label'_normd_amp/'$sta'_avgd_noisenormd_amp'
#head -2 $label'_normd_amp/'$sta'_avgd_normd_amp'
end
