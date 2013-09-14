cp /home/tianye/code/Script/GMT/gmtdefaults4 ./.gmtdefaults4
gmtset LABEL_OFFSET 0
gmtset HEADER_OFFSET 0

set fsrc = source.lst
set fsta = station_syn.lst
set outf = 'Source_Station.ps'
set REG = -R-40/40/-30/30
set SCA = -JX15
set title = 'Thetao2 Far Amp Sources'
psbasemap $REG $SCA -G255 -Ba10f1/a10f1:."$title":WSEn -K -P -X3 -Y7 >! $outf
set norm = `awk 'BEGIN{a=0}{if(a<$3)a=$3}END{print a}' $fsrc`
awk -v norm=$norm '{print $1,$2,$3/norm}' $fsrc | psxy -R -J -Sc -W2,red -O -K >> $outf
awk '{print $2,$3}' $fsta | psxy -R -J -St.3 -W1 -Gblue -O -K >> $outf
awk 'NR==1{print $2,$3, 10, 0, 1, 3, $1}' $fsta | pstext -R -J -O -K >> $outf
awk 'NR==6{print $2,$3, 10, 0, 1, 9, $1}' $fsta | pstext -R -J -O -K >> $outf

pwd | psxy -R -J -O >> $outf
echo $outf

