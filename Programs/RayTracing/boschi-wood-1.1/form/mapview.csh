#! /bin/csh
# arguments 
# 1 = name of xyz file to plot
# 2 = grid spacing
gmtset BASEMAP_TYPE plain LABEL_FONT_SIZE 12 LABEL_FONT 6 ANOT_FONT_SIZE 10 ANOT_FONT 6	
xyz2grd $1 -Gimage.bin -I$2/$2 -R0/360/-90/90 -L -N0.0 -V
grdimage image.bin -P -X0.5 -Y20 -R -Bg45/g45 -JH180/14 -K -V -Crel.cpt >! $1.ps
psscale -Crel.cpt -D15/3.5/4/.5 -O -L -K -B0.03:"relative anomaly": >> $1.ps
pscoast -R -O -JH -W1 -Dc -A10000 >> $1.ps
alias rm rm -r
rm image.bin
