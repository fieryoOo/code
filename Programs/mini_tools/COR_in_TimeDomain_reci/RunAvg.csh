#!/bin/csh
if ( $#argv != 2 )then
echo "usage: "$argv[0]" [input_sac] [output_sac]"
exit
endif

if( ! -e $argv[1] ) then
echo "File "$argv[1]" not exist!"
exit
endif

sac << END
r $argv[1]
abs
smooth mean h 128
smooth mean h 128
w temp.sm
r $argv[1]
divf temp.sm
w $argv[2]
quit
END

rm -f temp.sm
