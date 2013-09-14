#/bin/csh
if ($#argv != 2)then
echo "usage: amp_diff [sac_path_1] [sta_loc_2]"
exit
endif

cd $argv[1]
#rm -f event_TA_XR.lst
#foreach event (`ls -d 20*/ | cut -d/ -f1`)
#ls -d $argv[2]/20* | awk -F/ '{print $NF}' | grep $event >> event_TA_XR.lst
#end

cp /home/tianye/code/Script/GMT/gmtdefaults4 ./.gmtdefaults4

foreach station (`awk '{if(NR>1){print}}' $argv[2] | sed s/' '/'@'/g | sed s/'\+'/''/g`)
set sta_name=`echo $station | cut -d@ -f1`
set long=`echo $station | cut -d@ -f2`
set lati=`echo $station | cut -d@ -f3`
if (`echo "$long*2000000/2" | bc`<0)then
 set long=`echo "$long+360" | bc`
endif
#set gplong=`echo "$long*200/2*0.01" | bc`
#set gplati=`echo "$lati*200/2*0.01" | bc`
set longlow=`echo "$long*20/2*0.1" | bc | sed s/'\.0'/''/g`
set longhigh=`echo "$long*20/2*0.1+0.1" | bc | sed s/'\.0'/''/g`
set latilow=`echo "$lati*20/2*0.1" | bc | sed s/'\.0'/''/g`
if(`echo "$latilow*20/2" | bc` < 0)set latilow= `echo "$latilow-0.1" | bc | sed s/'\.0'/''/g`
set latihigh=`echo "$latilow+0.1" | bc | sed s/'\.0'/''/g`
#echo $long,$longlow,$longhigh;echo $lati,$latilow,$latihigh
set dis_lowlow=`echo $longlow $long $latilow $lati | awk '{printf "%.10f\n",($1-$2)^2+($3-$4)^2+0.0000000001}'`
set dis_lowhigh=`echo $longlow $long $latihigh $lati | awk '{printf "%.10f\n",($1-$2)^2+($3-$4)^2+0.0000000001}'`
set dis_highlow=`echo $longhigh $long $latilow $lati | awk '{printf "%.10f\n",($1-$2)^2+($3-$4)^2+0.0000000001}'`
set dis_highhigh=`echo $longhigh $long $latihigh $lati | awk '{printf "%.10f\n",($1-$2)^2+($3-$4)^2+0.0000000001}'`
foreach period (22 30 40 50 60 70 80 90 100)
echo 'working on station '$sta_name' for '$period'sec...'
set dis=`echo "$period*12" | bc`
cd $period'sec_10snr_'$dis'dis'
mkdir -p $argv[1]/$period'sec_10snr_'$dis'dis'/amp_diff_$period's_v2'
rm -f $argv[1]/$period'sec_10snr_'$dis'dis'/amp_diff_$period's_v2'/$sta_name
rm -f $argv[1]/$period'sec_10snr_'$dis'dis'/amp_diff_$period's_v2'/$sta_name'_pctg'
foreach event (`more ../event_TA_XR.lst`)
#set event=20080819163018
if ( ! -e $event'_am.txt_v2' ) then
continue
endif
if ( ! `awk -v long=$long -v lati=$lati '(long-$1)^2+(lati-$2)^2<0.00001' $argv[1]/$period'sec_10snr_'$dis'dis'/$event'_am.txt_v2' | wc -l` )then
continue
endif

awk -v long=$long -v lati=$lati 'sqrt(1.0000*(($1-long)^2+($2-lati)^2))>0.44 && sqrt(1.0000*(($1-long)^2+($2-lati)^2))<2.26' $event'_am.txt_v2' > throw_v2_onesta_am_temp
if ( `more throw_v2_onesta_am_temp | wc -l` < 10 )continue

set rg_longlow=`echo "$long*20/2*0.1-2.2" | bc`
set rg_longhigh=`echo "$long*20/2*0.1+2.3" | bc`
set rg_latilow=`echo "$lati*20/2*0.1-2.2" | bc`
set rg_latihigh=`echo "$lati*20/2*0.1+2.3" | bc`
echo '-R'$rg_longlow'/'$rg_longhigh'/'$rg_latilow'/'$rg_latihigh > temp_onesta_am.region
#echo $file'_no_'$sta_name
/home/tianye/code/Script/GMT/C_plot_travel throw_v2_onesta_am_temp temp_onesta_am.region

set amp_sta2=`awk -v long=$long -v lati=$lati '(long-$1)^2+(lati-$2)^2<0.00001 {print $3}' $argv[1]/$period'sec_10snr_'$dis'dis'/$event'_am.txt_v2'`
set lowlow=`grep $longlow'	'$latilow'	' throw_v2_onesta_am_temp.HD | awk '{print $3}'`
set lowhigh=`grep $longlow'	'$latihigh'	' throw_v2_onesta_am_temp.HD | awk '{print $3}'`
set highlow=`grep $longhigh'	'$latilow'	' throw_v2_onesta_am_temp.HD | awk '{print $3}'`
set highhigh=`grep $longhigh'	'$latihigh'	' throw_v2_onesta_am_temp.HD | awk '{print $3}'`
#echo $lowlow,$lowhigh,$highlow,$highhigh

set amp_pdt=`echo $dis_lowlow $lowlow $dis_lowhigh $lowhigh $dis_highlow $highlow $dis_highhigh $highhigh | awk '{printf "%.5f\n",(1/$1*$2+1/$3*$4+1/$5*$6+1/$7*$8)/(1/$1+1/$3+1/$5+1/$7)}'`
#echo $gplong,$gplati
#echo $amp_sta2 $amp_pdt
set amp_diff=`echo "$amp_sta2-($amp_pdt)" | bc`
set diff_pctg=`echo "$amp_diff*1000000/$amp_pdt*0.0001" | bc`

echo $amp_diff >> $argv[1]/$period'sec_10snr_'$dis'dis'/amp_diff_$period's_v2'/$sta_name
echo $diff_pctg >> $argv[1]/$period'sec_10snr_'$dis'dis'/amp_diff_$period's_v2'/$sta_name'_pctg'
#exit
end
#exit
cd ..
end
end
