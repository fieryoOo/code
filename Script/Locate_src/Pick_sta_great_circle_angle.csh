#!/bin/csh
if ($#argv != 6)then
echo "USAGE: Pick_sta_line.csh [long1] [lati1] [azimuth] [max_width] [station.lst] [out_staion.lst]"
exit
endif

awk '{print $2,$3,$1}' $argv[5] | project -C$argv[1]'/'$argv[2] -A$argv[3] -Q -S -W-$argv[4]'/'$argv[4] -Fzxy > $argv[6]
