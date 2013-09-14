#!/bin/csh
cp ~/.gmtdefaults4 .
gmtset HEADER_FONT = 4
gmtset HEADER_FONT_SIZE = 14
gmtset HEADER_OFFSET = 0.07c
gmtset ANNOT_FONT_PRIMARY = 4
gmtset ANNOT_FONT_SECONDARY = 4
gmtset LABEL_FONT = 4
gmtset ANNOT_OFFSET_PRIMARY = 0.08c
gmtset ANNOT_OFFSET_SECONDARY = 0.08c
gmtset LABEL_FONT_SIZE = 12p
gmtset ANNOT_FONT_SIZE_PRIMARY = 14p
gmtset LABEL_OFFSET = 0.12c

gmtset PLOT_DEGREE_FORMAT = +ddd:mm:ss
gmtset BASEMAP_TYPE = plain

set fcpt = vel.cpt
set Qmod = Qg
set REG = -R0.5/3.5/11/60
set SCA = -JX9/-5.2
set xlabel = "Age (Ma)"
set ylabel = "Depth (km)"
set outps = 'Models_Obs_HSC_MELT_NOMELT_'${Qmod}'.ps'

#plot observed model
rm -f observed.mod
foreach age ( 0.5 1.0 1.5 2.0 2.5 3.5 )
set fmod = MC.1.${age}.mod
awk '{print '$age',$1,$2}' 'Observed/'$fmod >> observed.mod
end
set title = "Observed"
awk '{if ($3>3.9) print $1,$2,$3}' observed.mod | surface -T0.2 $REG -G1.grd -I0.02/0.1
psbasemap -Ba0.5f0.1:."$title"::"$xlabel":/a10f5:"$ylabel":WSne -X4 -Y14 $SCA -R -K >! $outps
grdimage 1.grd -C$fcpt -R -J -O -K >> $outps
grdcontour 1.grd -C0.1 -W2 -R -J -O -K >> $outps
grdcontour 1.grd -C0.05 -W1,0/0/0,- -R -J -O -K >> $outps

#plot Half-space-cooling model
rm -f HSC.mod
foreach age ( 0.5 1.0 1.5 2.0 2.5 3.0 3.5 )
set fmod = CC.${age}
awk '{print '$age',$1,$2}' 'HSC/'$fmod >> HSC.mod
end
set title = "HalfSpaceCooling"
awk '{if ($3>3.9) print $1,$2,$3}' HSC.mod | surface -T0.2 -R -G1.grd -I0.02/0.1
psbasemap -Ba0.5f0.1:."$title"::"$xlabel":/a10f5:"$ylabel":WSne -X12 $SCA -R -O -K >> $outps
grdimage 1.grd -C$fcpt -R -J -O -K >> $outps
grdcontour 1.grd -C0.1 -W2 -R -J -O -K >> $outps
grdcontour 1.grd -C0.05 -W1,0/0/0,- -R -J -O -K >> $outps

#plot Goes model with melt
#awk '{print $4,$2,$3}' separate_${Qmod}Melt1s/Dis_Dep_Vs_Age.txt > Melt.mod
awk '$1<80 && $3<4{print $3,$1,$5}' separate_${Qmod}Melt1s/Seis_DehydDepl${Qmod}Melt1s.out > Melt.mod
#~/code/Programs/Smoothing/Gauss_Smoothing_map Melt.mod 0.5
set title = "${Qmod}Melt"
surface Melt.mod -T0.2 -R -G1.grd -I0.004/0.1
#grdsample 1.grd -Q -G2.grd -R -I0.02/0.2
psbasemap -Ba0.5f0.1:."$title"::"$xlabel":/a10f5:"$ylabel":WSne -X-12 -Y-8 $SCA -R -O -K >> $outps
grdimage 1.grd -C$fcpt -R -J -O -K >> $outps
grdcontour 1.grd -C0.1 -W2 -R -J -O -K >> $outps
grdcontour 1.grd -C0.05 -W1,0/0/0,- -R -J -O -K >> $outps

#plot Goes model with melt
awk '$1<80 && $3<4{print $3,$1,$5}' separate_${Qmod}1s/Seis_DehydDepl${Qmod}1s.out > NoMelt.mod
set title = "${Qmod}NoMelt"
surface NoMelt.mod -T0.2 -R -G1.grd -I0.004/0.1
psbasemap -Ba0.5f0.1:."$title"::"$xlabel":/a10f5:"$ylabel":WSne -X12 $SCA -R -O -K >> $outps
grdimage 1.grd -C$fcpt -R -J -O -K >> $outps
grdcontour 1.grd -C0.1 -W2 -R -J -O -K >> $outps
grdcontour 1.grd -C0.05 -W1,0/0/0,- -R -J -O -K >> $outps

#plot color bar
psscale  -C$fcpt -D-2/-2/10/0.3h -O -K -L  >> $outps

pwd | psxy -R -J -O >> $outps
echo $outps
