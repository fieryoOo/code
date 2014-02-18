#! /bin/csh
# arguments 
# 1 = name of xyz file with x component
# 2 = name of xyz file with y component
# 3 = grid spacing
gmtset HEADER_FONT_SIZE 20 OBLIQUE_ANOTATION 0 DEGREE_FORMAT 0 BASEMAP_TYPE plain
xyz2grd $1 -Gx.grd -I$3/$3 -R0/360/-90/90 -L -N0.0 -V
xyz2grd $2 -Gy.grd -I$3/$3 -R0/360/-90/90 -L -N0.0 -V
grdvector x.grd y.grd -E -K -P -X0.5 -Y20 -Bg45/g45 -JQ180/14 -V -W2 -S.05 >! $1.ps 
psxy /home/boschil/Documents/GMT/plate_boundaries.codes -W3/150/150/150 -M -O -K -V -R -JQ >> $1.ps
pscoast -R -O -JQ -W1 -Dc -A10000 >> $1.ps
alias rm rm -r
rm x.grd y.grd
