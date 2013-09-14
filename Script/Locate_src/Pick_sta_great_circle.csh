#!/bin/csh
if ($#argv != 7)then
echo "USAGE: Pick_sta_line.csh [long1] [lati1] [long2] [lati2] [max_width] [station.lst] [out_staion.lst]"
exit
endif

awk '{print $2,$3,$1}' $argv[6] | project -C$argv[1]'/'$argv[2] -E$argv[3]'/'$argv[4] -Lw -Q -S -W-$argv[5]'/'$argv[5] -Fzxy > $argv[7]
