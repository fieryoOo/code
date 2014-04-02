#!/bin/csh
if ($#argv != 2) then
  echo "USAGE: "$0" [in_seedf_list] [label]"
  exit 1
endif

awk -F$argv[2]_ '{print $2,$1}' $argv[1] | awk -F. '{print $1,$2,$3,$0}' | sed s/'JAN'/'1'/ | sed s/'FEB'/'2'/ | sed s/'MAR'/'3'/ | sed s/'APR'/'4'/ | sed s/'MAY'/'5'/ | sed s/'JUN'/'6'/ | sed s/'JUL'/'7'/ | sed s/'AUG'/'8'/ | sed s/'SEP'/'9'/ | sed s/'OCT'/'10'/ | sed s/'NOV'/'11'/ | sed s/'DEC'/'12'/ | awk -v lab=$argv[2] '{print $5lab"_"$4,$1,$2,$3}' | sort -s -g -k4 | sort -s -g -k3 | sort -s -g -k2
