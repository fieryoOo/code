#!/bin/csh
if ( $#argv != 3 )then
echo "usage: Avg_l_normd_Amp.csh [sac_path] [sta_p_dist_days_table] [sta_lst]"
exit
endif

cd $argv[1]
set label=`echo $argv[3] | awk -F/ '{print $NF}' | cut -d. -f1 | cut -d_ -f2`
#echo "touching those fucking foolish SAC files..."
awk '{print $1}' sta_for_amp_norm_days.lst > sta_for_amp_norm.lst
#touch_sac sta_for_amp_norm.lst
#echo "Normalizing the amplitude for each station pairs..."
#/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/amp_length_normalize sta_for_amp_norm_days.lst

mkdir -p $label'_bin_normd_amp'
foreach sta_info (`awk '{print $1"@"$2"@"$3}' $argv[3]`)
set sta=`echo $sta_info | cut -d@ -f1`
set long=`echo $sta_info | cut -d@ -f2`
set lati=`echo $sta_info | cut -d@ -f3`
rm -f $sta'_file.lst'
foreach file ( `grep '_'$sta sta_for_amp_norm.lst | awk '{print $1"_amp_l_norm.txt"}'` )
if(-e $file)echo $file >> $sta'_file.lst'
end
if ( ! -e $sta'_file.lst' || `more $sta'_file.lst' | wc -l` < 8 )then
echo $sta": less than 8 stations!"
rm -f $sta'_file.lst'
rm -f $label'_bin_normd_amp/'$sta'_'*'_daynormd_amp'
rm -f $label'_bin_normd_amp/'$sta'_'*'_noisenormd_amp'
continue
endif
echo "Extracting the normd amplitude...Working on station "$sta"..."
set per_file=`head -1 $sta'_file.lst'`
awk '{print $1}' $per_file > $label'_bin_normd_amp/per_file'
foreach period (`more $label'_bin_normd_amp/per_file'`)
echo $sta $long $lati > $label'_bin_normd_amp/'$sta'_'$period'_daynormd_amp'
echo $sta $long $lati > $label'_bin_normd_amp/'$sta'_'$period'_noisenormd_amp'
foreach file (`more $sta'_file.lst'`)
set sta2=`echo $file | sed s/'_'$sta/''/ | cut -d_ -f2 | cut -d. -f1`
set daynorm_amp=`awk -v sta2=$sta2 -v per=$period '$1==per {print $2}' $file`
set noisenorm_amp=`awk -v sta2=$sta2 -v per=$period '$1==per {print $3}' $file`
echo $sta2 $daynorm_amp >> $label'_bin_normd_amp/'$sta'_'$period'_daynormd_amp'
echo $sta2 $noisenorm_amp >> $label'_bin_normd_amp/'$sta'_'$period'_noisenormd_amp'
end
end
end

cd $argv[1]'/'$label'_bin_normd_amp'
foreach period (`more per_file`)
ls *$period'_daynormd_amp' > 'daynormd_'$period'.lst'
ls *$period'_noisenormd_amp' > 'noisenormd_'$period'.lst'
#/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/bin_avg_center_sta 'daynormd_'$period'.lst' /home/tianye/data_days/station_all.lst
#/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/bin_avg_center_sta 'noisenormd_'$period'.lst' /home/tianye/data_days/station_all.lst
/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/azi_cov_center_sta 'daynormd_'$period'.lst' /home/tianye/data_days/station_all.lst
/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/azi_cov_center_sta 'noisenormd_'$period'.lst' /home/tianye/data_days/station_all.lst
end

cd $argv[1]'/'$label'_bin_normd_amp'
foreach ctg (day noise)
foreach period (`more per_file`)
rm -f $ctg'normd_'$period'.txt'
echo "Predicting the normd amp for each sta...Working on "$ctg" "$period"sec..."
foreach file (`ls *$period'_'$ctg'normd_amp.txt'`)
/home/tianye/code/Script/GMT/C_plot_travel $file /home/tianye/data_days/region_TA
set long=`head -1 $file | awk '{print $1}'`
set lati=`head -1 $file | awk '{print $2}'`
set amp=`awk -v long=$long -v lati=$lati '$1<long+0.1 && $1>long-0.1 && $2>lati-0.1 && $2<lati+0.1' $file'.HD' | awk 'BEGIN{a=0;weight=0}{w=1/sqrt((long-$1)**2+(lati-$2)**2);a+=$3*w;weight+=w}END{print a/weight}'`
echo $long $lati $amp >> $ctg'normd_'$period'.txt'
end
end
end
