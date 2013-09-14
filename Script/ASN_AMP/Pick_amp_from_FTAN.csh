#!/bin/csh
if ( $#argv != 2 )then
echo "usage: Pick_amp_from_FTAN.csh [COR_path] [sta_lst]"
exit
endif

cd $argv[1]

foreach per (7 15)
set i = 0
foreach sta1 (`awk '{print $1}' $argv[2]`)
@ i += 1
cd $sta1
foreach sta2 (`awk -v i=$i 'NR>i {print $1}' $argv[2]`)
set amp=`awk -v per=$per 'BEGIN{p=-99999;a=0}{if(per>p && per<$3){amp=(a/(per-p)+$6/($3-per))/(1/(per-p)+1/($3-per))};p=$3;a=$6}END{print amp}' 'COR_'$sta1'_'$sta2'.SAC_2_DISP.1'`
echo 'per: '$per'  sta: '$sta1'-'$sta2'  amp: '$amp
end
cd ..
end
end
