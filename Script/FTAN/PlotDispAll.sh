#!/bin/bash

if [ $# != 1 ]; then
	echo "Usage: "$0" [disp file list]"
	exit
fi

# ps file
psout=${1}.ps
pwd | psxy -R -J -K > $psout

# plot group
psbasemap -R2/35/0.5/4.5 -JX8 -B:."Group":WeSn -X2 -Y3 -O -K >> $psout
while read fdisp stmp; do
	awk '{print $3,$4}' $fdisp | psxy -R -J -Sc0.1 -Gred -O -K >> $psout
done < $1

# plot phase
psbasemap -R2/35/0.5/4.5 -JX8 -B:."Phase":WeSn -X12 -O -K >> $psout
while read fdisp stmp; do
	awk '{print $3,$5}' $fdisp | psxy -R -J -Sc0.1 -Gsteelblue -O -K >> $psout
done < $1

# end
pwd | psxy -R -J -O >> $psout
