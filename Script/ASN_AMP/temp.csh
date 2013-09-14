cd /utera/tianye/data_avg_monthly_stack
foreach sta (LRL VES X25A)
echo "Combining data for station "$sta"..."
set dir1=`ls -d */ | grep -m1 $sta | cut -d/ -f1`
set file1=`ls $dir1'/'$dir1'_stacked_'*'_days.SAC_amp_l_norm.txt' | head -1`
foreach period (`awk '{print $1}' $file1`)
rm -f combined_$sta'_'$period'_noisesqr_days_noise_ampm_ampc.txt'
foreach pair (`ls -d */ | grep $sta`)
more $pair'/'$period'_noisesqr_days_noise_ampm_ampc.txt' >> combined_$sta'_'$period'_noisesqr_days_noise_ampm_ampc.txt'
end
end
end


foreach file ( `ls -d combined_*_noisesqr_days_noise_ampm_ampc.txt` )
set file_new=`echo $file | sed s/'ampc'/'ampc_2nd_year'/g`
awk '$2>280' $file > $file_new
end
