#!/bin/csh
#echo on
if ($#argv != 2) then
  echo "USAGE: calc_sta_p_time_length.csh [sta_p_lst] [rec_path]"
  exit 1
endif

cd $argv[2]
set out_name_1=`echo $argv[1] | awk -F/ '{print $NF}' | cut -d. -f1`
set out_name=`echo ${out_name_1}"_withdays.txt"`
rm -f ${out_name}
foreach pair (`awk '{print $1"@"$2"@"$3}' $argv[1]`)
set sta1=`echo $pair | cut -d@ -f1`
set sta2=`echo $pair | cut -d@ -f2`
set dist=`echo $pair | cut -d@ -f3`
set dir=$sta1'_'$sta2
echo "working on station pair "$dir
set day_num=0
foreach year (2008 2009 2010)
@ day_num += `/home/tianye/code/Programs/CORRECT_ASN_AMP/calc_cor_time rec_$year'/'$sta1'_rec.SAC' rec_$year'/'$sta2'_rec.SAC' 2`
end
echo $sta1 $sta2 $dist $day_num >> ${out_name}
end

