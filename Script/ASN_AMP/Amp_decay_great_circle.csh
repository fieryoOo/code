#!/bin/csh
if ( $#argv != 6 )then
echo "usage: Amp_decay_great_circle.csh [COR_dir (station.lst in it)] [source_long] [source_lati] [azimuth] [max_width] [per]"
exit
endif

cd $argv[1]
#set per=16.070479
#set per=10.889532
set per=$argv[6]
set fstep=-0.08
set outf=$argv[2]'_'$argv[3]'_'$argv[4]
awk '{print $2,$3,$1}' station.lst | project -C$argv[2]'/'$argv[3] -A$argv[4] -Q -S -W-$argv[5]'/'$argv[5] -Fzxy > $outf'_sta.lst'
if( `more $outf'_sta.lst' | wc -l` < 2 )then
echo "less than 2 station found!"
rm -f $outf'_sta.lst'
exit
endif

foreach info (`awk '{print $1"@"NR"@"$2"@"$3}' $outf'_sta.lst'`)
set sta1=`echo $info | cut -d@ -f1`
set staN=`echo $info | cut -d@ -f2`
set long1=`echo $info | cut -d@ -f3`
set lati1=`echo $info | cut -d@ -f4`
rm -f $outf'_'$sta1'_dist_amp'
grep $sta1 file_daynum.lst > temp_daynum.lst
foreach info2 (`awk -v staN=$staN 'NR>staN {print $1"@"$2"@"$3}' $outf'_sta.lst'`)
set sta2=`echo $info2 | cut -d@ -f1`
#set file=$sta1'/COR_'$sta1'_'$sta2'.SAC'
set file=$sta1'/COR_'$sta1'_'$sta2'.SAC_amp_snr'
echo "file: "$file
if( -e $file ) then
set filec=$sta1'/COR_'$sta1'_'$sta2'.SAC'
if( `awk -v file=$filec '$1==file {if($2<100){print 1}}' temp_daynum.lst` )continue
awk -v per=$per -v step=$fstep 'BEGIN{f=1/per;fl=1/exp(log(f)-step);fh=1/exp(log(f)+step)} $1<fh && $1>fl && ( ($3>5 && $5>8) || $3>8) {print ($1-per)**2,$2,$3}' $file | sort -g | head -2 > amp_temp
else
set file=$sta2'/COR_'$sta2'_'$sta1'.SAC_amp_snr'
if( ! -e $file )then
echo "no info between "$sta1" and "$sta2". Skipped!"
continue
endif
set filec=$sta2'/COR_'$sta2'_'$sta1'.SAC'
if( `awk -v file=$filec '$1==file {if($2<100){print 1}}' temp_daynum.lst` )continue
awk -v per=$per -v step=$fstep 'BEGIN{f=1/per;fl=1/exp(log(f)-step);fh=1/exp(log(f)+step)} $1<fh && $1>fl && ( ($5>5 && $3>8) || $5>8) {print ($1-per)**2,$4,$5}' $file | sort -g | head -2 > amp_temp
endif
if( `more amp_temp | wc -l` < 2 ) continue
set amp=`awk 'BEGIN{a=0;w=0}{a=a+$2/($1+1e-10);w=w+1/($1+1e-10)}END{print a/w}' amp_temp`
set snr=`awk 'BEGIN{a=0;w=0}{a=a+$3/($1+1e-10);w=w+1/($1+1e-10)}END{print a/w}' amp_temp`
#if( ! -e $file ) then
#set file2=$sta2'/COR_'$sta2'_'$sta1'.SAC'
#echo "file2: "$file2
#if( ! -e $file2 ) then
#echo "no info between "$sta1" and "$sta2". Skipped!"
#continue
#endif
#sac << END
#r $file2
#reverse
#w $file
#quit
#END
#endif
set long2=`echo $info2 | cut -d@ -f2`
set lati2=`echo $info2 | cut -d@ -f3`
set dist=`/home/tianye/code/Programs/DIST/get_dist $lati1 $long1 $lati2 $long2`
echo $dist $amp $snr $sta2 >> $outf'_'$sta1'_dist_amp'
end #foreach info2
if (`more $outf'_'$sta1'_dist_amp' | wc -l` < 2)then
rm -f $outf'_'$sta1'_dist_amp'
endif

#awk '{print $1}' $outf'_'$sta1'.lst' > temp.lst
#foreach fname (`more temp.lst`)
#rm -f $fname'_amp_l_norm_pos.txt'
#end
#/home/tianye/code/Programs/CORRECT_ASN_AMP_DAYNUM/amp_length_normalize temp.lst pos
#rm -f temp.lst $outf'_'$sta1'_dist_amp'
#foreach infof (`awk '{print $1"@"$2}' $outf'_'$sta1'.lst'`)
#set file=`echo $infof | cut -d@ -f1`'_amp_l_norm_pos.txt'
#if( ! -e $file )continue
#set sta2=`echo $infof | cut -d_ -f3 | cut -d. -f1`
#set dist=`echo $infof | cut -d@ -f2`
#set amp=`awk -v per=$per '($1-per)**2<0.05 {print $2}' $file`
#echo $dist $amp $sta2 >> $outf'_'$sta1'_dist_amp'
#end #foreach infof
end #foreach info
rm -rf $outf
mkdir $outf
mv $outf'_'* $outf
cd $outf
foreach ampf (`ls $outf'_'*'_dist_amp'`)
#if( `more $ampf | wc -l` < 2 )then
#rm -f $ampf
#else 
set csta=`echo $ampf | cut -d_ -f4`
grep $csta $outf'_sta.lst' >> temps.lst
#endif
end
awk '{print $3,$1,$2}' temps.lst | sort | awk '{print $2,$3,$1}' > $outf'_center_sta.lst'
rm -r temps.lst
cp '../USglobal_R_'$per'.HD'* .
/home/tianye/code/Script/GMT/C_plot_input_region_sta_1deg USglobal_R_$per.HD.cpt USglobal_R_$per.HD /home/tianye/region_TA $outf'_center_sta.lst'
#awk '{print "\"250_10_10_"$1"_dist_amp\" w l ,\\"}' $outf'_center_sta.lst'
awk -v outf=$outf '{printf "\""outf"_"$1"_dist_amp\""; if(NR<=5){printf " lt 1"} else if(NR<=10){printf " lt 2"} else if(NR<=15){printf " lt 3"} else if(NR<=20){printf " lt 4"} else {printf " lt 5"} printf ",\\\n"}' $outf'_center_sta.lst' > plot.txt
