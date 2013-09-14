if ($#argv != 1)then
echo "usage: Statistic_avg [sac_path]"
exit
endif

cd $argv[1]
#set table_file='/home/tianye/data_Eikonal/SAC_XR/amp_diff_avg_table_XR'
foreach period (30 40 50 60 70 80 90 100)
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis/amp_diff_'$period's/'
rm -f amp_diff_table_pctg_$period's'
#echo '\n'$period'sec' >> $table_file
foreach file (`ls SC??_pctg.picked`)
if (`more $file | wc -l` >= 10)then
#echo $file >> $table_file
set sta=`echo $file | cut -d_ -f1`
awk -v sta=$sta 'BEGIN{a=0}{a=a+$1}END{print sta,a/NR}' $file >> amp_diff_table_pctg_$period's'
endif
end
cd ../..
end
