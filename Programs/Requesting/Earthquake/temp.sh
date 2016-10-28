#!/bin/bash

while read sta lon lat; do
	awk -v sta=$sta '$1==sta' station.loc~ | head -n1
done < /work1/tianye/EQKLocation/station_avail.lst
