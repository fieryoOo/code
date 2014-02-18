#! /bin/csh
# arguments 
# 1 = name of xyz file with x component
# 2 = name of xyz file with y component
# 3 = name of xyz file with x component
# 4 = name of xyz file with y component
# 5 = grid spacing
gmtset HEADER_FONT_SIZE 20 OBLIQUE_ANOTATION 0 DEGREE_FORMAT 0 BASEMAP_TYPE plain
xyz2grd $1 -Gx1.grd -I$5/$5 -R0/360/-90/90 -L -N0.0 -V
xyz2grd $2 -Gy1.grd -I$5/$5 -R0/360/-90/90 -L -N0.0 -V
xyz2grd $3 -Gx2.grd -I$5/$5 -R0/360/-90/90 -L -N0.0 -V
xyz2grd $4 -Gy2.grd -I$5/$5 -R0/360/-90/90 -L -N0.0 -V
grdvector x1.grd y1.grd -E -K -P -X0.5 -Y20 -Bg45/g45 -JQ180/14 -V -W2 -S.02 >! $1.ps 
grdvector x2.grd y2.grd -E -O -K -P -Bg45/g45 -JQ180/14 -V -W2 -S.02 >> $1.ps 
pscoast -R -O -JQ -W1 -Dc -A10000 >> $1.ps
alias rm rm -r
rm x1.grd y1.grd x2.grd y2.grd

