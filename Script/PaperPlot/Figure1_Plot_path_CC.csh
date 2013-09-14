#!/bin/csh
# Script for plotting the inversion results for the western US region from the weekly stacks.

# GET PLOTTING PARAMETERS

### set gmt parameters
cp /home/tianye/code/Script/GMT/gmtdefaults4 ./.gmtdefaults4
#gmtset HEADER_FONT_SIZE 20
#gmtset ANNOT_FONT_SIZE 12
#gmtset HEADER_OFFSET -0.1
#gmtset PLOT_DEGREE_FORMAT -ddd:mm:ssF
#gmtset LABEL_OFFSET 0
#gmtset LABEL_FONT 4
#gmtset ANNOT_FONT 4
#gmtset HEADER_FONT 4

#set input_map = area/OBS_topo.HD
set cptfile = age/age.3.6.xyz.cpt
set input_map = age/age.3.6.xyz
set res = 0.01
set region_f = area/region_OBS
set sta_lst=`echo $input_map | awk -F.HD '{print $1}'`
set sta_name1 = 'area/OBS_topo_name1'
set sta_name2 = 'area/OBS_topo_name2'
set fpath = 'path.lst'

set title='Path Coverage'

set out_name = `echo $input_map `
set output_ps_file = ${out_name}.ps
echo "writing $output_ps_file file"
if (-f $output_ps_file) rm $output_ps_file

set REG = `more $region_f`
set long1=`echo $REG | cut -d/ -f1 | sed s/'\-R'/''/` 
set long2=`echo $REG | cut -d/ -f2`
set lat1=`echo $REG | cut -d/ -f3`
set lat2=`echo $REG | cut -d/ -f4`
set long=`echo $long1 $long2 | awk '{print ($1+$2)/2}'`

set latimin=`echo $REG | cut -d/ -f3`
set latimax=`echo $REG | cut -d/ -f4`
#set SCA = -JX5i/5i
set SCA = -Jm1.3

#pwd | psxy -H $REG $SCA -X3 -Y8.0  -V -K  >! $output_ps_file
#psbasemap $REG $SCA -X5 -Y6 -V -K >! $output_ps_file

#xyz2grd $input_map -Gtomo.grd -I$res -V $REG
#grdsample tomo.grd -Q -Gtomo2.grd $REG -I$sms
echo $input_map
surface $input_map $REG -T0.2 -Gtomo.grd -I0.01
#grdgradient tomo.grd -A0 -N2 -Gtomo_gradient.grd
### shift 1 ###
grdimage $SCA $REG -Ba2f0.5g1/a1f0.5g1WSen:."$title": tomo.grd -C$cptfile -X5 -Y6 -K >! $output_ps_file #-Itomo_gradient.grd

#pscoast $SCA $REG -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 -S135/206/235 >> $output_ps_file
#pscoast $SCA $REG -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 >> $output_ps_file

#psxy /home/tianye/code/Programs/head/wus_province_II.dat $SCA $REG -W5/255/0/0 -M"99999 99999"  -O -K >> $output_ps_file
psxy /home/tianye/code/Programs/head/platebound.gmt $SCA $REG -W8/70/120/160 -M"99999 99999"  -O -K >> $output_ps_file
psscale -B1.5:"Lithospheric Age (Ma)": -C$cptfile -D3.9/-1/7.8/0.3h -O -K  >> $output_ps_file

psxy $fpath $REG $SCA -W5/70/70/70 -m -O -K >> $output_ps_file
psxy $sta_name2 $SCA $REG -St.3 -W3 -Gwhite  -O -K >> $output_ps_file
#pstext $sta_name2 $SCA $REG -Gblack -O -K >> $output_ps_file
#awk '{printf "%f %f ",$1,$2}' $sta_name | psxy $REG $SCA -Svs0.03c/0.4c/0.15c -W4,red -Gred -O -K >> $output_ps_file
psxy $sta_name1 $REG $SCA -W6/255/0/60 -O -K >> $output_ps_file #0/80/255
psxy $sta_name1 $SCA $REG -St.3 -W3 -G255/0/60  -O -K >> $output_ps_file
pstext $sta_name1'_1' $SCA $REG -Gblack -Wwhite,o2,red -O -K >> $output_ps_file
#echo 230.0 46.06 | psxy $SCA $REG -Sa.3 -Gorangered -O -K >> $output_ps_file
#echo 229.9 46.1 12 0 1 3 CHS | pstext $SCA $REG -Gblack -O -K >> $output_ps_file
pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $output_ps_file
-0.3 8.3 20 0.0 4 LT (a)
EOF

rm -f tomo.grd tomo2.grd tomo_gradient.grd

###file names###
set sacf=CC/J29A_J47A/COR_J29A_J47A.SAC
set stap=`echo $sacf | cut -d/ -f1`
set disf1=$sacf'_1_DISP.0'
set disf2=$sacf'_2_DISP.1'
set ampf1=$sacf'_1_AMP'
set ampf2=$sacf'_2_AMP'
set snrf=$sacf'_amp_snr'
#set outf=$stap'_CC_FTAN.ps'
###check file existences###
foreach file ( $sacf $disf1 $disf2 $ampf1 $ampf2 $snrf )
if( ! -e $file )then
echo "File "$file" not found."
exit
endif
end
###plot time series###
set t_info=`saclst KEVNM KSTNM DIST KCMPNM f $sacf | awk '{print "ev/csta: "$2"  sta: "$3"  dist: "$4"  ch: "$5}'`
set vmin=0.5
set vmax=4.5
set tb=`saclst b f $sacf | awk '{print $2}'`
set te=`saclst e f $sacf | awk '{print $2}'`
set sacff=$sacf
set sym_flag=`echo $tb $te | awk '{printf "%.0f\n",$1/$2}'`
if ( $sym_flag == -1 ) then
set sacff=$sacf'_sym'
set tlen=`echo $tb $te | awk '{a=-$1;if(a>$2){a=$2} print a}'`
echo "tlen" $tlen
/home/tianye/Software/sac/bin/sac << END
cut -$tlen 0 
r $sacf
reverse
ch b 0
ch nzjday 1
w temp_n
cut 0 $tlen
r $sacf
ch nzjday 1
addf temp_n
w $sacff
quit
END
rm -f temp_n
endif
set sptf=$sacff'.am'
set dist=`saclst DIST f $sacf | awk '{print $2}'`
set wb=`echo $dist $vmin $vmax | awk '{printf "%.0f",$1*(1./$3-(1./$2-1./$3)/4.)}'`
if( $wb < 0 ) set wb=0
set we=`echo $dist $vmin $vmax | awk '{printf "%.0f",$1*(1./$2+(1./$2-1./$3)/4.)}'`
if( $we > $te ) set we=$te
### change sym or both
set wb = `echo $we | awk '{print -$1}'`
###
#set yrange=`echo $wb $we | awk '{printf "%.2g/%.2g",$1,$2}'`
#set ytic=`echo $wb $we | awk '{printf "%.1g",($2-$1)/20.}'`
#set ymrk=`echo $ytic | awk '{print $1*5}'`
set yrange='-500/500'
set ytic = 20
set ymrk = 200
set xtmp=`saclst DEPMIN DEPMAX f $sacff | awk '{a=-$2;if(a<$3){a=$3} print a}'`
#set xrange=`echo $xtmp | awk '{a=$1*1.1; printf "%.2g/%.2g",-a,a}'`
#set xtic=`echo $xtmp | awk '{printf "%.1g",$1/5.}'`
#set xmrk=`echo $xtmp | awk '{printf "%.2g",$1*1.1}'`
set xrange = '-0.00003/0.00003'
set xtic = 0.00001
#set xmrk = 0.00003
set xmrk = 0.00005
echo "x: "$xrange "  y: "$yrange

set REG='-R'$yrange'/'$xrange
### shift 2
#-Ba3f3:"Period (sec)":/a0.5f0.1:"Velocity (km/s)":/a0.5f0.1:."FTAN Diagram":WSne
   ### original sac
#/home/tianye/code/Programs/pssac/pssac $REG -JX3.5i/.8i -B'a'$ymrk'f'$ytic':'"Time (sec)"':/a'$xmrk'f'$xtic':'"Amplitude (nm/sec)"::."Time Series":'WeSn' -W1 -X11 -Y-0.5 $sacf -K -O >> $output_ps_file
/home/tianye/code/Programs/pssac/pssac $REG -JX3.75i/.85i -B'a'$ymrk'f'$ytic':'"Time (sec)"':/a'$xmrk'f'$xtic':'::."Time Series":'WeSn' -W1 -X10 -Y7.6 $sacf -K -O >> $output_ps_file
   ### negative 0.9 ~ 4.2
#set b = `echo $dist | awk '{print -$1/0.9}'`
#set e = `echo $dist | awk '{print -$1/4.2}'`
set b = 55
set e = 250
echo "distance: "$dist
echo "time window: "$b" "$e
echo $dist | awk '{print "vel range: "$1/55., $1/250.}'
sac << END
cut $b $e
r $sacf
w cut1
quit
END
   ### positive 0.9 ~ 4.2
#set b = `echo $dist | awk '{print $1/4.2}'`
#set e = `echo $dist | awk '{print $1/0.9}'`
set b = -250
set e = -55
sac << END
cut $b $e
r $sacf
w cut2
quit
END
/home/tianye/code/Programs/pssac/pssac $REG -JX3.75i/.85i cut1 -W,red -K -O >> $output_ps_file
/home/tianye/code/Programs/pssac/pssac $REG -JX3.75i/.85i cut2 -W0/80/255 -K -O >> $output_ps_file
pstext -R0/10/0/10 -JX3.5i/.8i -V -O -K -N << EOF >>  $output_ps_file
-0.5 16.2 20 0.0 4 LT (b)
EOF
rm -f cut1 cut2
#/home/tianye/code/Programs/pssac/pssac $REG -JX3.5i/.8i -B'a'$ymrk'f'$ytic':'"Time (s)"':/a'$xmrk'f'$xtic':'"Amplitude (nm/sec)"::."Time Series":'WeSn' -X10 $sacf -K -O >> $output_ps_file
#/home/tianye/code/Programs/pssac/pssac  -JX3.5i/.8i -B'a'$ymrk'f'$ytic':'"Time (s)"':/a'$xmrk'f'$xtic':'"Amplitude (nm/sec)"::."Time Series":'WeSn' -X10 $sacf -K -O >> $output_ps_file
#pstext -R0/10/0/10 -JX6.5i/1.5i -O -K -N << END >> $output_ps_file
#-1.7 4.6 12 270 20 LT time(s)
#12 0 15 270 20 LT $t_info
#END

###plot amplitude and group dispersion 2###
set perl=`awk 'BEGIN{a=999}{if(a>$3)a=$3}END{print a}' $disf2`
set perh=`awk 'BEGIN{a=0}{if(a<$3)a=$3}END{print a}' $disf2`
set REG=`echo $vmin $vmax $perl $perh | awk '{print "-R"log($3)"/"log($4)"/"$1"/"$2}'`
echo $REG > region_tmp
set SCA = -JX3.5i/2.2i
awk '{print log($1),$2,$3}' $ampf2 > amp_tmp
#TXT2CPT_FTAN amp_tmp
C_plot_travel_nn amp_tmp region_tmp 0.02 0.2
xyz2grd amp_tmp.HD -Gamp.grd -I0.02 -V $REG
### shift 3
grdimage $REG $SCA amp.grd -X0.5 -Y-8.4 -C'CC/amp.cpt' -O -K >> $output_ps_file
set REG="-R"$perl"/"$perh"/"$vmin"/"$vmax
set SCAl = `echo $SCA | sed s/'i'/'il'/`
#rm -f amp_tmp amp_tmp.cpt amp_tmp.HD
#awk '{print $4,$3}' $disf1 | psxy $REG -JX2.5i/-2.8il -Ba0.5f0.1:"group velocity (km/sec)":/a3f3:."Dispersion(2nd)":WSne -Sc.15 -W6,white -A -O -K >> $output_ps_file
awk '{print $3,$4}' $disf2 | psxy $REG $SCAl -Ba3f3:"Period (sec)":/a1.f0.2+0.5:"Velocity (km/sec)":/a0.5f0.1:."FTAN Diagram":WSne -Sc.12 -Gred -W4 -A -O -K >> $output_ps_file
awk '{print $3,$5}' $disf2 | psxy $REG $SCAl -Sc.12 -W4,black -Gwhite -A -O -K >> $output_ps_file
pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $output_ps_file
-0.5 12.2 20 0.0 4 LT (c)
EOF

#echo -1.5 2.5 20 0.0 4 LT Figure 2 | pstext -R0/1/0/1 $SCA -V -O -K -N >>  $output_ps_file
pwd | psxy -R -J -O >> $output_ps_file
echo $output_ps_file
