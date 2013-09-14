#!/bin/csh 
# Script to run the resolution analysis version of tomography
limit datasize 1000000; limit stacksize 120000

# GET PERIOD FOR RUNNING TOMOGRAPHY
if ($#argv != 1) then
  echo "USAGE: $0 [period]"
  exit 1
endif
set per=$argv[1]
#echo "Running itomo_ra_sp_cu_shn at ${per} seconds"
echo "Running itomo_ra_shn_l at ${per} seconds"
set time1resid=`date | tr -s ' ' | cut -d' ' -f2,3,4`

# CHECK FOR EXISTENCE OF data[period]s_ph.txt FILE
set data="data/data${per}s_ph.txt"
if ( !(-f $data) ) then
  echo "$data file does not exist"
  exit 1
endif

if (!(-d $per)) mkdir $per
if (-f contour.ctr)  rm contour.ctr
#cp /data/bassoon/morganm/TOMO/reqd_files/contour.ctr .
cp /home/morganm/tomo_files/contour.ctr .

set alpha=400
#set alpha=300
# beta value (aka, alpha-1)
set beta=100
set sigma=70
set name="USArr_"$alpha"_"$sigma

#/home/morganm/bin/itomo_ra_sp_cu_shn $data $name $per << EOF
/home/morganm/bin/itomo_ra_sp_cu_shn_l $data $name $per << EOF
me
4
5 
$beta
6
$alpha
$sigma
$sigma
7
20 55 0.5
8
230 250 0.5
12
.05
0.5
16
19
v
q
go
EOF

set dir = $per"/"$alpha"_"$sigma
if(!(-d $dir)) mkdir $dir
mv USArr* $dir

# SEND EMAIL TO NOTIFY OF COMPLETION
set time2resid=`date | tr -s ' ' | cut -d' ' -f2,3,4`
mail -s "$0 - $per" morganm@ciei.colorado.edu << EOF
$0 complete: $per second result
`pwd`

start time: $time1resid
end time: $time2resid

EOF

exit
