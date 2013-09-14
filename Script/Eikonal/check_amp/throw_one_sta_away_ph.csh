#!/bin/csh
if ( $#argv != 2 )then
echo "usage: throw_sta.csh [sac_path] [sta_list_for_throw]"
exit
endif

cd $argv[1]
foreach station (`awk '{if(NR>0){print}}' $argv[2] | sed s/' '/'@'/g | sed s/'\+'/''/g`)
set sta_name=`echo $station | cut -d@ -f1`
set long=`echo $station | cut -d@ -f2`
set lati=`echo $station | cut -d@ -f3`
if (`echo "$long*2/2" | bc`<0)then
 set long=`echo "$long+360" | bc`
endif
set gplong=`echo "$long*200/2*0.01" | bc`
set gplati=`echo "$lati*200/2*0.01" | bc`

foreach period (22 30 40 50 60 70 80 90 100)
echo "working on "$sta_name" for "$period"sec..."
set dis=`echo "$period*12" | bc`
cd $period'sec_10snr_'$dis'dis'
rm -rf ph_v2_onesta_threw_$sta_name
mkdir -p ph_v2_onesta_threw_$sta_name

foreach file (`ls 20*.ph.txt_v2`)
#set file=20090418191758.ph.txt
if ( `awk -v long=$long -v lati=$lati '$1<long+2 && $1>long-2 && $2<lati+2 && $2>lati-2' $file | wc -l` < 10 )continue
if (`grep $gplong $file | grep -c $gplati`)then
set gptemp=`grep $gplong $file | grep $gplati`
awk -v long=$long -v lati=$lati '$1<long+2 && $1>long-2 && $2<lati+2 && $2>lati-2' $file | grep -v "$gptemp" > ph_v2_onesta_threw_$sta_name'/'$file'_no_'$sta_name
else
awk -v long=$long -v lati=$lati '$1<long+2 && $1>long-2 && $2<lati+2 && $2>lati-2' $file > ph_v2_onesta_threw_$sta_name'/'$file'_no_'$sta_name
endif
set longlow=`echo "$long*20/2*0.1-2" | bc`
set longhigh=`echo "$long*20/2*0.1+2.1" | bc`
set latilow=`echo "$lati*20/2*0.1-2" | bc`
set latihigh=`echo "$lati*20/2*0.1+2.1" | bc`
echo '-R'$longlow'/'$longhigh'/'$latilow'/'$latihigh > temp.region
#echo $file'_no_'$sta_name
/home/tianye/code/Script/GMT/C_plot_travel ph_v2_onesta_threw_$sta_name'/'$file'_no_'$sta_name temp.region
end
cd ..
end

end
