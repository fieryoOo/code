#!/bin/bash
lstdir=./seed.lst
year=(2001 2004 2005 2006 2007 2008 2009)
month=(1 2 3 4 5 6 7 8 9 10 11 12)
day=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)
if [ -e seed.lst ]; then
        rm -f seed.lst
fi
for yy in ${year[@]}; do
	if [ $[$yy % 400] -eq "0" ]; then
		leapflag=1
      	elif [ $[$yy % 4] -eq 0 ]; then
		  if [ $[$yy % 100] -ne 0 ]; then
			  leapflag=1
      		  else
			  leapflag=0
    		  fi
	else
		leapflag=0
	fi
	for mm in ${month[@]}; do
		for dd in ${day[@]}; do
			if [ $dd -ge 30 ] && [ $mm -eq 2 ]; then
				continue
			elif [ $dd -eq 31 ] && ( [ $mm -eq 4 ] || [ $mm -eq 6 ] || [ $mm -eq 9 ] || [ $mm -eq 11 ] ) ; then
				continue
			elif [ $mm -eq 2 ] && [ $dd -eq 29 ] && [ $leapflag -eq 0 ]; then
				continue
			fi
				echo ${yy}_${mm}_${dd}.seed $yy $mm $dd >> seed.lst
		done
        done
done
mv seed.lst $lstdir
