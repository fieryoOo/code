#!/bin/csh
if ($#argv != 1) then
  echo "USAGE: Plot_FTAN [sac_file_name]"
  exit 1
endif

set symf = $argv[1]'_2_DISP.1'
set posf = $argv[1]'_pos_2_DISP.1'
set negf = $argv[1]'_neg_2_DISP.1'
foreach file ($symf $posf $negf)
if(! -e $file) then
echo "Cannot open file "$file". stopped."
exit
endif
end
set outf = $argv[1]'_ph.ps'
gnuplot << END
set term postscript enhanced color
set out '$outf'
plot '$symf' using 3:5 lt -1, '$posf' using 3:5 lt 1, '$negf' using 3:5 lt 3
set term x11
set out
q
END
echo $outf
