#!/bin/csh
#echo on
if ($#argv != 1) then
  echo "USAGE: tool_do_AtoD  period"
  exit 1
endif

set per=$argv[1]

@ dis = 12 * $per
cd $per"sec_10snr_"$dis"dis"
ls 20*ph.txt_v2 | awk '{print "~/Programs/C_plot_travel_am",$1}' > plot_all_am.csh
csh plot_all_am.csh
/home/tianye/code/Programs/amp_HD_to_amp_gradient_HD_to_amp_laplace_HD event_lst $per
cd ..

