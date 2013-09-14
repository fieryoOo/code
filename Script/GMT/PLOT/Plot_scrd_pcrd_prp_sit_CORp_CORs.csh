#!/bin/csh
if ($#argv != 2) then
  echo "USAGE: Plot_scrd_pcrd_prp_sit_CORp_CORs.csh [per] [0(non-smth) or 1(smthd)]"
  exit 1
endif

set per=$argv[1]
set sm_h=70
cp /home/tianye/code/Script/GMT/gmtdefaults4 ./.gmtdefaults4
gmtset HEADER_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 8
gmtset HEADER_OFFSET 0.1

set cpt_f = ( '/media/WORK/tianye/Events_Attenuation/alpha_cpt/alpha_'$per'.cpt' '/media/WORK/tianye/Events_Attenuation/alpha_cpt/alpha_'$per'_rvs.cpt' '/media/WORK/tianye/Events_Attenuation/alpha_cpt/alpha_'$per'_pos.cpt' )
if( $argv[2] ) then
set f1='alpha_sit_crd_map_1_smd_'$sm_h'.HD'
set f2='alpha_prp_crd_map_1_smd_'$sm_h'.HD'
set f3='alpha_prp_map_1_smd_'$sm_h'.HD'
set f4='alpha_site_map_1_smd_'$sm_h'.HD'
set output_ps_file = 'scrd_prp_pcrd_sit_CORp_CORs_smd_'$sm_h'.ps'
else
set f1='alpha_sit_crd_map_1.HD'
set f2='alpha_prp_crd_map_1.HD'
set f3='alpha_prp_map_1.HD'
set f4='alpha_site_map_1.HD'
set output_ps_file = 'scrd_prp_pcrd_sit_CORp_CORs.ps'
endif

set REG=`more ~/region_TA`
set long1=`echo $REG | cut -d/ -f1 | sed s/'\-R'/''/`
set long2=`echo $REG | cut -d/ -f2`
set long=`echo $long1 $long2 | awk '{print ($1+$2)/2}'`
set latimin=`echo $REG | cut -d/ -f3`
set latimax=`echo $REG | cut -d/ -f4`
set SCA = -JN$long/2.5i
set res=0.1

if (-f $output_ps_file) rm $output_ps_file
rm -f temp.lst

set avg=( 0 0 0 0 )
set nf = 1
foreach input_map ( $f1 $f2 $f3 $f4 )

set sta_lst=`echo $input_map | awk -F.HD '{print $1}'`
awk '{print $1"_"$2,$3}' $sta_lst > $sta_lst'_tmp'
echo $sta_lst'_tmp' >> temp.lst
set avg[$nf]=`awk 'BEGIN{a=0}{a+=$3}END{print a/NR}' $sta_lst`
awk -v avg=$avg[$nf] '{print $1,$2,$3-avg}' $input_map > temp1

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
xyz2grd temp1 -Gtomo.grd -I$res -V $REG
grdsample tomo.grd -Q -Gtomo2.grd $REG -I$sms
grdimage $SCA $REG tomo2.grd -C$cptfile -Ba5f1/a3f1WESN:."$title": -O -P -K >> $output_ps_file

pscoast $SCA $REG -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 -S135/206/235 >> $output_ps_file

psxy /home/tianye/code/Programs/head/wus_province_II.dat $SCA $REG -W5/255/0/0 -M"99999 99999"  -O -K >> $output_ps_file

psxy $sta_lst $SCA $REG -St.1 -W3 -Ggray -O -K >> $output_ps_file

@ nf ++
end #foreach

psscale  -C$cpt_f[1] -P -D-12.3/7/10/0.5 -O -K -L  >> $output_ps_file

# LABELS
#pstext -R0/10/0/10 -JX10c  -V -O -N -G0 -Y3c -X0.5c << END >>  $output_ps_file
#3 -4.5 16 0.0 7 CB  event: $event_label
#END

\rm tomo.grd tomo2.grd topo_westUS2.xyz

~/code/Programs/mini_tools/File_Combiner temp.lst 1 temp1
set outf='loc_raw_prp_sit_pcrd_scrd_crd.txt'
if( $argv[2] ) set outf='loc_raw_prp_sit_pcrd_scrd_crd_smd_'$sm_h'.txt'
awk '{print $1,$2-$8,$6,$8,$4,$2,$2+$6}' temp1 > $outf
rm -f *_tmp temp.lst temp1

if( $per == 40 ) then
set REG=-R-0.0015/0.0015/-0.0015/0.0015
else if ( $per == 60 ) then
set REG=-R-0.001/0.001/-0.001/0.001
else
set REG=-R-0.002/0.002/-0.002/0.002
endif

set SCA=-JX2.5i/2.5i
gmtset HEADER_OFFSET -0.5

awk -v avg3=$avg[3] -v avg1=$avg[1] '{print $3-avg3, $6-avg1}' $outf > temp1
set slp1=`/home/tianye/code/Programs/FIT/least_squares_line_origin temp1 0 | awk '{printf "%.3f",$1}'`
set slp2=`/home/tianye/code/Programs/FIT/least_squares_line_origin temp1 1 | awk '{printf "%.3f",$1}'`
set slp=`echo $slp1 $slp2 | awk '{print 0.5*(atan2($1,1)+atan2($2,1))}' | awk '{printf "%.3f",sin($1)/cos($1)}'`
set corr=`/home/tianye/code/Programs/mini_tools/Correlation_infile temp1 | awk '{printf "%.3f",$1}'`
pwd | psxy -H $REG $SCA -X-11 -Y-9.5 -P -V -O -K  >> $output_ps_file
psxy temp1 -B0.001/0.001:."scrd v.s. prp": $REG $SCA -Sc.05 -Gred -O -K >> $output_ps_file
echo $slp | awk '{print "-1",-$1,"\n1",$1}' | psxy $REG $SCA -Gblack -A -O -K >> $output_ps_file
pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $output_ps_file
5.1 9.5 10 0.0 4 LT  slp: $slp1 ~ $slp2
7 8.8 10 0.0 4 LT  cslp: $slp
7 8.1 10 0.0 4 LT  corr: $corr
EOF
awk -v avg4=$avg[4] -v avg2=$avg[2] '{print $4-avg4, $5-avg2}' $outf > temp1
set slp1=`/home/tianye/code/Programs/FIT/least_squares_line_origin temp1 0 | awk '{printf "%.3f",$1}'`
set slp2=`/home/tianye/code/Programs/FIT/least_squares_line_origin temp1 1 | awk '{printf "%.3f",$1}'`
set slp=`echo $slp1 $slp2 | awk '{print 0.5*(atan2($1,1)+atan2($2,1))}' | awk '{printf "%.3f",sin($1)/cos($1)}'`
set corr=`/home/tianye/code/Programs/mini_tools/Correlation_infile temp1 | awk '{printf "%.3f",$1}'`
pwd | psxy -H $REG $SCA -X10 -P -V -O -K  >> $output_ps_file
psxy temp1 -B0.001/0.001:."pcrd v.s. sit": $REG $SCA -Sc.05 -Gred -O -K >> $output_ps_file
echo $slp | awk '{print "-1",-$1,"\n1",$1}' | psxy $REG $SCA -Gblack -A -O -K >> $output_ps_file
pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $output_ps_file
5.1 9.5 10 0.0 4 LT  slp: $slp1 ~ $slp2
7 8.8 10 0.0 4 LT  cslp: $slp
7 8.1 10 0.0 4 LT  corr: $corr
EOF

rm -f temp1

