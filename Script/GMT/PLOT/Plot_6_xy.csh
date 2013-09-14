#!/bin/csh
if ($#argv != 2) then
  echo "USAGE: Plot_raw_crd_prp_sit_HGr_HGc.csh [in_file.lst (f1 f2) (max_row: 6)] [outf_name]"
  exit 1
endif

cp /home/tianye/code/Script/GMT/gmtdefaults4 ./.gmtdefaults4
gmtset HEADER_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 6.5
gmtset HEADER_OFFSET -0.5

set REG='-R5/30/0.5/4.5'
set SCA = -JX3.3il/3i

set output_ps_file = $argv[2]
if (-f $output_ps_file) rm $output_ps_file

set Nf = `more $argv[1] | wc -l`
set nf = 1
foreach input ( `awk '{print $1"@"$2}' $argv[1]` )
set f1=`echo $input | cut -d@ -f1`
set f2=`echo $input | cut -d@ -f2`

set title='node '$nf

if( $nf % 2 == 0 ) then
set X0=10
set Y0=0
else
set X0=-10
set Y0=-9
endif

if( $nf == 1 )then
pwd | psxy -H $REG $SCA -X11.5 -Y30 -P  -V -K  >! $output_ps_file
else pwd | psxy -H $REG $SCA -X$X0 -Y$Y0 -P  -V -O -K  >> $output_ps_file
endif

echo "5 3.8\n30 3.8" | psxy $REG $SCA -W1,black -O -K >> $output_ps_file
echo "5 3.5\n30 3.5" | psxy $REG $SCA -W1,black -O -K >> $output_ps_file
psxy $f1 -Ba3f3/0.5:."$title":SnWe $REG $SCA -W3,red -O -K >> $output_ps_file
psxy $f2 -Ba3f3/0.5:."$title":SnWe $REG $SCA -W3,blue -O -K >> $output_ps_file

@ nf ++
end #foreach


