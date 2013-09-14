if ( $#argv != 2 )then
echo "usage: throw_sta.csh [sac_path] [sta_list_for_throw] [region_infile]"
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

foreach period (22)
echo "working on "$sta_name" for "$period"sec..."
set dis=`echo "$period*12" | bc`
cd $period'sec_10snr_'$dis'dis'
mkdir am_v2_onesta_threw_$sta_name

foreach file (`ls 20*_am.txt_v2`)
#set file=20090418191758_am.txt_v2
if (`grep $gplong $file | grep -c $gplati`)then
grep -v "`grep $gplong $file | grep $gplati`" $file > am_v2_onesta_threw_$sta_name'/'$file'_no_'$sta_name
else
cp $file am_v2_onesta_threw_$sta_name'/'$file'_no_'$sta_name
endif
/home/tianye/code/Script/GMT/C_plot_travel am_v2_onesta_threw_$sta_name'/'$file'_no_'$sta_name argv[3]
end
cd ..
end

end
