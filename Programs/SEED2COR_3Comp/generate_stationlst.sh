#!/bin/bash
seeddir=/rc_scratch/life9360/BOSA_COR/SEED
#seeddir=/lustre/janus_scratch/life9360/COR_US_CONTINENT/SEED
#lstdir=/lustre/janus_scratch/life9360/COR_US_CONTINENT/working_2004
lstdir=/lustre/janus_scratch/life9360/COR_BOSA
if [ -e station.lst ]; then
	        rm -f station.lst
fi
cd $seeddir
rm s1.lst
ls *seed >> seed_fname.lst
for sname in `cat seed_fname.lst`
do
	rdseed -S -f $sname
        cat rdseed.stations >> s1.lst
done
sort -k1,1 -u s1.lst > s2.lst
awk '{print $1,$4,$3}' s2.lst > station.lst
rm s1.lst s2.lst seed_fname.lst
cp station.lst $lstdir
