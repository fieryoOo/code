#!/bin/csh
if ($#argv != 3) then
  echo "USAGE: Plot_raw_crd_prp_sit_HGr_HGc.csh [in_file.lst (HD cpt) (max: 6)] [outf_name] [region_infile]"
  exit 1
endif

cp /home/tianye/code/Script/GMT/gmtdefaults4 ./.gmtdefaults4
gmtset HEADER_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 6.5
gmtset HEADER_OFFSET 0.1

set REG=`more $argv[3]`
set long1=`echo $REG | cut -d/ -f1 | sed s/'\-R'/''/`
set long2=`echo $REG | cut -d/ -f2`
set long=`echo $long1 $long2 | awk '{print ($1+$2)/2}'`
set latimin=`echo $REG | cut -d/ -f3`
set latimax=`echo $REG | cut -d/ -f4`
set SCA = -JN$long/2.5i
set res=0.1

set output_ps_file = $argv[2]
if (-f $output_ps_file) rm $output_ps_file

set Nf = `more $argv[1] | wc -l`
set nf = 1
foreach input ( `awk '{print $1"@"$2}' $argv[1]` )
set input_map=`echo $input | cut -d@ -f1`
set cptfile=`echo $input | cut -d@ -f2`
set sta_lst=`echo $input_map | awk -F.HD '{print $1}'`

foreach i (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
set event_label=`echo $input_map | sed s/'\.'/'_'/g | cut -d_ -f$i`
if( `echo $event_label | awk '{print length}'` == 14 )then
break
endif
end

set title=`echo $input_map | awk -F/ '{print $NF}' | sed s/'\.'/'_'/g | sed s/'HD'/''/g | sed s/'_'/' '/g`

if( $nf % 2 == 0 ) then
set X0=10
set Y0=0
else
set X0=-10
set Y0=-9
endif

if( $nf == 1 )then
pwd | psxy -H $REG $SCA -X12 -Y30.7 -P  -V -K  >! $output_ps_file
else pwd | psxy -H $REG $SCA -X$X0 -Y$Y0 -P  -V -O -K  >> $output_ps_file
endif

xyz2grd $input_map -Gtomo.grd -I$res -V $REG
set sms=`echo $res | awk '{print $1/5.0}'`
grdsample tomo.grd -Q -Gtomo2.grd $REG -I$sms
grdimage $SCA $REG tomo2.grd -C$cptfile -Ba5f1/a3f1WESN:."$title": -O -P -K >> $output_ps_file

if( $nf == 1 )then                                            <
set cont_dis=$per                                             <
grdcontour $SCA $REG tomo2.grd -C$cont_dis -W.5p,black -O -K  <
endif

pscoast $SCA $REG -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 -S135/206/235 >> $output_ps_file

psxy /home/tianye/code/Programs/head/wus_province_II.dat $SCA $REG -W5/255/0/0 -M"99999 99999"  -O -K >> $output_ps_file

psxy $sta_lst $SCA $REG -St.1 -W3 -Ggray -O -K >> $output_ps_file
if( $nf == $Nf ) then
psscale  -C$cptfile -E -I -D3.2/-1/6.4/0.15h -O -L  >> $output_ps_file
else
psscale  -C$cptfile  -E -I -D3.2/-1/6.4/0.15h -O -K -L  >> $output_ps_file
endif

@ nf ++
end #foreach

\rm tomo.grd tomo2.grd topo_westUS2.xyz

