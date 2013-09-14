#extract misfit for each station
set infile = sta1_sta2_dist_pvelp_pveln_misfit.txt
set stalst = station_18.loc
foreach sta (`awk '{print $3}' $stalst`)
set outf = $sta'_misfit.txt'
awk -v sta=$sta '$1==sta{print $6,$2,$3,$4,$5}' $infile > $outf
awk -v sta=$sta '$2==sta{print -$6,$1,$3,$5,$4}' $infile >> $outf
if( `more $outf | wc -l` == 0 ) rm -f $outf
end #sta

cp /home/tianye/code/Script/GMT/gmtdefaults4 ./.gmtdefaults4
#compute statistics and plot histograms
#gmtset HEADER_FONT_SIZE 30
#gmtset ANNOT_FONT_SIZE 20
#gmtset LABEL_FONT_SIZE 25
#gmtset HEADER_OFFSET 0
#gmtset PLOT_DEGREE_FORMAT -ddd:mm:ssF
#gmtset LABEL_OFFSET 0
#gmtset ANNOT_OFFSET 0.2
#gmtset LABEL_FONT 4
#gmtset ANNOT_FONT 4
#gmtset HEADER_FONT 4

rm -f statistics.txt
foreach info (`awk '{print $3"_misfit.txt@"$1"@"$2}' $stalst`)
set file = `echo $info | cut -d@ -f1`
set lon = `echo $info | cut -d@ -f2`
set lat = `echo $info | cut -d@ -f3`
set mean = `awk 'BEGIN{a=0}{a+=$1}END{print a/NR}' $file`
set std = `awk 'BEGIN{std=0}{std+=($1-'$mean')**2}END{print (std/(NR-1+1e-3))**0.5}' $file`
set N = `more $file | wc -l`
echo $file $mean $std $N $lon $lat >> statistics.txt
pshistogram -JX6/10 -R-6/6/0/6 -Ba5f1:"time (sec)":/a5f1:"# of measurements"::."mean $mean  std $std  # $N":WSEn -P -X3 -Y3 -W1 -F -Z0 -L3,red $file > $file'.ps'
end

set outf = statistics.txt.ps
### plot histogram of all measurement misfits
set mean = `awk 'BEGIN{a=0}{a+=$6}END{printf "%.2g",a/NR}' $infile`
set stdd = `awk 'BEGIN{std=0}{std+=($6-'$mean')**2}END{printf "%.5g",(std/(NR-1+1e-3))**0.5}' $infile`
set std = `echo $stdd | awk '{printf "%.2g", $1}'`
set N = `more $infile | wc -l`
awk '{print $6}' $infile | pshistogram -JX6/9 -R-5/5/0/38 -Ba5f1:"time (sec)":/a5f1:"No. of Inter-station Paths"::."Phase Time Misfit":WSen -X7 -Y5 -W1 -F -Z0 -L3,red -K > $outf
echo "1 9.5 15 0.0 4 LT mean "$mean"   std "$std | pstext -R0/10/0/10 -J -V -O -K -N  >>  $outf # "  # "$N
echo "-1.8 11 20 0.0 4 LT (a)" | pstext -R0/10/0/10 -J -V -O -K -N  >>  $outf

### plot mean misfit on each station (dots)
#awk 'NR>1{print $2}' statistics.txt | pshistogram -JX6/9 -R-2/2/0/10.5 -Ba1f0.5:"time (sec)":/a2f1:."Station Mean Misfits":WSen -X8 -Y-1 -W0.4 -F -Z0 -L3,red -K -O >> $outf
awk 'NR>1{print $2,$4}' statistics.txt | sort -g | psxy -Ba0.5f0.1:"time (sec)":/a2f1:"No. of Measurements Per Station"::."Station Mean Misfits":WSen -JX6/9 -R-1./1./2/15 -Sc.15 -W4/0/80/250 -X9  -O -K >> $outf

### plot 1-sigma confidence interval
#echo "-0.77 0\n-0.77 15\n>\n0.77 0\n0.77 15" | psxy -R -J -W3,red -m -O -K >> $outf
set N = 1
rm -f plot.tmp
while( $N < 32 )
echo $stdd $N | awk '{print $1/sqrt($2/2.), $2/2.}' >> plot.tmp
@ N ++
end
awk '{print $1,$2}' plot.tmp | psxy -R -J -W2,red,- -O -K >> $outf
awk '{print -$1,$2}' plot.tmp | psxy -R -J -W2,red,- -O -K >> $outf
awk '{print $1*2,$2}' plot.tmp | psxy -R -J -W2,red -O -K >> $outf
awk '{print -$1*2,$2}' plot.tmp | psxy -R -J -W2,red -O -K >> $outf
pslegend -R0/10/0/5 -D7/4.45/2.3/1.2/BL -J -O -K -F1 -G255 << EOF >> $outf
S .2i - .3i - 1p,red,4:0p .5i  1-@~s
S .2i - .3i - 1p,red .5i  2-@~s
EOF

echo "-2.2 11.1 20 0.0 4 LT (b)" | pstext -R0/10/0/10 -J -V -O -K -N  >>  $outf
#awk 'NR>1{print $2,$4}' statistics.txt | sort -g | psxy -Ba0.5f0.1g0.77:"time (sec)":/a2f1:"Number of Measurements"::."Station Mean Misfits":WSen -JX6/9 -R-1/1/2/16 -Sc.1 -W2/0/80/250 -X9  -O -K >> $outf
# or plot mean misfit on each station (map)
if( 0 ) then
set cptfile = area/OBS_topo.cpt
set input_map = area/OBS_topo.HD
set res = 0.01
set region_f = area/region_OBS
set sta_lst=`echo $input_map | awk -F.HD '{print $1}'`
set title='Station Mean Misfits'
set REG = `more $region_f`
set SCA = -Jm1.1
set sms=`echo $res | awk '{print $1/5.0}'`
xyz2grd $input_map -Gtomo.grd -I$res -V $REG
grdgradient tomo.grd -A0 -N2 -Gtomo_gradient.grd
grdimage $SCA $REG -Ba2f0.5g1/a1f0.5g1WSen:."$title": tomo.grd -Itomo_gradient.grd -C$cptfile -X9 -Y1 -K -O >> $outf
psxy /home/tianye/code/Programs/head/platebound.gmt $SCA $REG -W8/70/70/70 -M"99999 99999"  -O -K >> $outf
awk '{print $5,$6,$2}' statistics.txt | psxy -Cpolar.cpt $SCA $REG -Sc.3 -O -K >> $outf
psscale -B0.4:"Misfit (sec)": -Cpolar.cpt -D3/-1/6/0.3h -O -K  >> $outf
echo "-2.3 8 20 0.0 4 LT (b)" | pstext -R0/10/0/10 -J -V -O -K -N  >>  $outf
endif

echo "-5 10 20 0.0 4 LT Figure 3" | pstext -R0/1/0/1 -J -V -O -N  >>  $outf
echo $outf
