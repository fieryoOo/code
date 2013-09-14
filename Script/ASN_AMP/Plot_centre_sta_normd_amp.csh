#!/bin/csh
if ( $#argv != 2 )then
echo "usage: Plot_centre_sta_normd_amp.csh [centre_sta] [max_dist]"
exit
endif

set per=12
set fstep=-0.08
cd /media/WORK/tianye/US_ASN_amp/STACK_08to11_norm
#grep -m1 $argv[1] station.lst > $argv[1]'.lst'

echo "Producing normd amp map..."
rm -f $argv[1]'/'$argv[1]'_center_normd_amp_'$per'sec.txt'
foreach file ( `grep $argv[1]'\/' file_daynum.lst | awk '$2>60 {print $1"_amp_snr"}'` )
awk -v per=$per -v step=$fstep 'BEGIN{f=1/per;fl=1/exp(log(f)-step);fh=1/exp(log(f)+step)} $1<fh && $1>fl && ( ($3>3 && $5>5) || $3>5) {print ($1-per)**2,$2,$3}' $file | sort -g | head -2 > amp_temp
if( `more amp_temp | wc -l` < 2 ) continue
set amp=`awk 'BEGIN{a=0;w=0}{a=a+$2/($1+1e-10);w=w+1/($1+1e-10)}END{print a/w}' amp_temp`
set sta2=`echo $file | cut -d_ -f3 | cut -d. -f1`
set location=`grep -m1 $sta2 station.lst | awk '{long=$2;if(long<0){long=long+360} print long,$3}'`
echo $location $amp >> $argv[1]'/'$argv[1]'_center_normd_amp_'$per'sec.txt'
end # file
foreach file ( `grep '_'$argv[1]'.SAC' file_daynum.lst | awk '$2>60 {print $1"_amp_snr"}'` )
awk -v per=$per -v step=$fstep 'BEGIN{f=1/per;fl=1/exp(log(f)-step);fh=1/exp(log(f)+step)} $1<fh && $1>fl && ( ($5>3 && $3>5) || $5>5) {print ($1-per)**2,$4,$5}' $file | sort -g | head -2 > amp_temp
if( `more amp_temp | wc -l` < 2 ) continue
set amp=`awk 'BEGIN{a=0;w=0}{a=a+$2/($1+1e-10);w=w+1/($1+1e-10)}END{print a/w}' amp_temp`
set sta2=`echo $file | cut -d_ -f2`
set location=`grep -m1 $sta2 station.lst | awk '{long=$2;if(long<0){long=long+360} print long,$3}'`
echo $location $amp >> $argv[1]'/'$argv[1]'_center_normd_amp_'$per'sec.txt'
end # file
#foreach finfo (`grep $argv[1]'.SAC' file_daynum.lst | awk '$2>60 {print $1"@"$2}'`)
#set file=`echo $finfo | cut -d@ -f1`
#set sta2=`echo $file | cut -d/ -f1`
#set dnum=`echo $finfo | cut -d@ -f2`
#set ofile=$argv[1]'/COR_'$argv[1]'_'$sta2'.SAC'
#sac << END
#r $file
#REVERSE
#w $ofile
#quit
#END
#echo $ofile'   '$dnum >> sta_for_amp_norm_days.lst
#end #finfo
#awk '{print $1}' sta_for_amp_norm_days.lst > SAC.lst
#touch_sac SAC.lst

#echo "Normalizing the amplitude for each station pairs..."
#rm -f $argv[1]'/COR_'$argv[1]*'.SAC_amp_l_norm_'$argv[3]'.txt'
#/home/tianye/code/Programs/CORRECT_ASN_AMP_DAYNUM/amp_length_normalize sta_for_amp_norm_days.lst $argv[3]

cd $argv[1]
echo "Plotting..."
#cp ../station.lst $argv[1]'.lst'
#set per_file=`ls 'COR_'$argv[1]'_'*'.SAC_amp_l_norm_'$argv[3]'.txt' | head -1`
#foreach period (`awk '{print $1}' $per_file`)
#echo 'Extracting the normd_amp for each period. Working on '$period'sec'
#rm -f $argv[1]'_center_normd_amp_'$argv[3]'_'$period'sec.txt'
#foreach file (`ls 'COR_'$argv[1]'_'*'.SAC_amp_l_norm_'$argv[3]'.txt'`)
#set sta2=`echo $file | cut -d_ -f3 | cut -d. -f1`
#set location=`grep -m1 $sta2 $argv[1]'.lst' | awk '{long=$2;if(long<0){long=long+360} print long,$3}'`
#set normd_amp=`awk -v per=$period '$1==per {print $2}' $file`
#set snr=`awk -v per=$period '$1==per {print $5}' $file`
#echo $location $normd_amp $snr >> $argv[1]'_center_normd_amp_'$argv[3]'_'$period'sec.txt'
#end
#end

#minmax $argv[1]'.lst' | awk '{print "-R/"$6"/"$7}' | sed s/'<'/''/g | sed s/'>'/''/g > $argv[1]'.region'

grep -m1 $argv[1] ../station.lst | awk '{long=$2;if(long<0){long=long+360} print $1,long,$3}' > $argv[1]'.loc'
set long=`awk '{print $2}' $argv[1]'.loc'`
set lati=`awk '{print $3}' $argv[1]'.loc'`
#set per_file=`ls 'COR_'$argv[1]'_'*'.SAC_amp_l_norm_'$argv[3]'.txt' | head -1`
#foreach period (`awk '{print $1}' $per_file`)
$argv[1]'/'$argv[1]'_center_normd_amp_'$per'sec.txt'
/home/tianye/code/Script/GMT/C_plot_travel_positive $argv[1]'_center_normd_amp_'$per'sec.txt' ~/region_TA
awk -v long=$long -v lati=$lati '$1>long-8 && $1<long+8 && $2>lati-10 && $2<lati+10' $argv[1]'_center_normd_amp_'$per'sec.txt.HD' > temp_cpt
/home/tianye/code/Script/GMT/TXT2CPT temp_cpt
mv temp_cpt.cpt $argv[1]'_center_normd_amp_'$per'sec.txt.HD.cpt'
rm -f temp_cpt
/home/tianye/code/Script/GMT/C_plot_kernel_input_sta $argv[1]'_center_normd_amp_'$per'sec.txt.HD.cpt' $argv[1]'_center_normd_amp_'$per'sec.txt.HD' ~/region_TA $argv[1]'.loc'
#end
