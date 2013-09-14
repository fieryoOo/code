#!/bin/csh
if ($#argv != 1) then
  echo "USAGE: Plot_raw_crd_prp_sit_HGr_HGc.csh [per]"
  exit 1
endif

set per=$argv[1]
cp /home/tianye/code/Script/GMT/gmtdefaults4 ./.gmtdefaults4
gmtset HEADER_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 8
gmtset HEADER_OFFSET 0.1

set cpt_f = ( '/media/WORK/tianye/Events_Attenuation/alpha_cpt/alpha_'$per'.cpt' '/media/WORK/tianye/Events_Attenuation/alpha_cpt/alpha_'$per'_rvs.cpt' '/media/WORK/tianye/Events_Attenuation/alpha_cpt/alpha_'$per'_pos.cpt' )
set f1='alpha_raw_map_1.HD'
set f2='alpha_crd_map_1.HD'
set f3='alpha_prp_map_1.HD'
set f4='alpha_site_map_1.HD'
set REG=`more ~/region_TA`
set long1=`echo $REG | cut -d/ -f1 | sed s/'\-R'/''/`
set long2=`echo $REG | cut -d/ -f2`
set long=`echo $long1 $long2 | awk '{print ($1+$2)/2}'`
set latimin=`echo $REG | cut -d/ -f3`
set latimax=`echo $REG | cut -d/ -f4`
set SCA = -JN$long/2.5i
set res=0.1

set output_ps_file = 'raw_crd_prp_sit_HGp_HGs.ps'
if (-f $output_ps_file) rm $output_ps_file
rm -f temp.lst

set nf = 1
foreach input_map ( $f1 $f2 $f3 $f4 )

set sta_lst=`echo $input_map | awk -F.HD '{print $1}'`

foreach i (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
set event_label=`echo $input_map | sed s/'\.'/'_'/g | cut -d_ -f$i`
if( `echo $event_label | awk '{print length}'` == 14 )then
break
endif
end

set title=`echo $input_map | sed s/'\.'/'_'/g | cut -d_ --complement -f$i | sed s/'_txt.*HD'/''/g | sed s/'_'/' '/g`

if( $nf % 2 == 0 ) then
set X0=9
set Y0=0
else
set X0=-9
set Y0=-8
endif

if( $nf == 1 )then
pwd | psxy -H $REG $SCA -X13.2 -Y29 -P  -V -K  >! $output_ps_file
else pwd | psxy -H $REG $SCA -X$X0 -Y$Y0 -P  -V -O -K  >> $output_ps_file
endif

set ncpt=`echo $nf | awk '{print int(($1+1)/2)}'`
set cptfile = $cpt_f[$ncpt]
set sms=`echo $res | awk '{print $1/5.0}'`
xyz2grd $input_map -Gtomo.grd -I$res -V $REG
grdsample tomo.grd -Q -Gtomo2.grd $REG -I$sms
grdimage $SCA $REG tomo2.grd -C$cptfile -Ba5f1/a3f1WESN:."$title": -O -P -K >> $output_ps_file

pscoast $SCA $REG -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 -S135/206/235 >> $output_ps_file

psxy /home/tianye/code/Programs/head/wus_province_II.dat $SCA $REG -W5/255/0/0 -M"99999 99999"  -O -K >> $output_ps_file

#psxy $sta_lst $SCA $REG -St.1 -W3 -Ggray -O -K >> $output_ps_file

@ nf ++
end #foreach

psscale  -C$cpt_f[1] -P -D-12.3/7/10/0.5 -O -K -L  >> $output_ps_file

# LABELS
#pstext -R0/10/0/10 -JX10c  -V -O -N -G0 -Y3c -X0.5c << END >>  $output_ps_file
#3 -4.5 16 0.0 7 CB  event: $event_label
#END

\rm tomo.grd tomo2.grd topo_westUS2.xyz

if( $per == 40 ) then
set REG=-R-0.0018/0.002/0/35
set width=0.00012
else if ( $per == 60 ) then
set REG=-R-0.0009/0.001/0/35
set width=0.00008
else
set REG=-R-0.0027/0.003/0/35
set width=0.0002
endif

set SCA=-JX3i/2.5i
gmtset HEADER_OFFSET -0.5

set file=`echo $f1 | cut -d. -f1`
set avg=`awk 'BEGIN{a=0}{a+=$3}END{print a/NR}' $file`
set std=`awk -v avg=$avg 'BEGIN{std=0}{std+=($3-avg)**2}END{print sqrt(std/(NR-1))}' $file`
awk '{print $3}' $file | pshistogram -B0.001/5:."raw  mean $avg  std $std": $REG -Z1 $SCA -W$width -L1.0p,red -K -Y-9.5 -X-11.5 -O -K >> $output_ps_file
set file=`echo $f2 | cut -d. -f1`
set avg=`awk 'BEGIN{a=0}{a+=$3}END{print a/NR}' $file`
set std=`awk -v avg=$avg 'BEGIN{std=0}{std+=($3-avg)**2}END{print sqrt(std/(NR-1))}' $file`
awk '{print $3}' $file | pshistogram -B0.001/5:."crd  mean $avg  std $std": $REG -Z1 $SCA -W$width -L1.0p,red -K -X10 -O >> $output_ps_file

