#!/bin/csh
if ( $#argv != 2)then
echo "USAGE: get_event [SAC_path] [sta.lst]"
exit
endif
cd $argv[1]
foreach per ( 6 8 10 14 18 24 30)
#foreach per ( 8 10 14 20 30 )

@ dis = 6 * $per
/home/tianye/code/Programs/ASN/ASN_gr_ph_amp_map $argv[2] $argv[2] file_daynum.lst $per 100 $dis

end
