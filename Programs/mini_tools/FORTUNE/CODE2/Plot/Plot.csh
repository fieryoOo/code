set fS = 'S.txt'
set fE = 'edges.txt'
set fout = 'VoronoiDiagram.ps'

cd ../mine/save_final
make clean
make
cd ../../Plot
../mine/save_final/Fortune $fS > $fE
#awk '{print $3}' out.txt | awk -F"->" '{print $1"\n"$2"\n>"}' | sed s/'('/''/g | sed s/')'/''/g | sed s/','/' '/g > $fE

set REG = -R0/100/0/100
set SCA = -JX15
psbasemap $REG $SCA -Ba10f5:"x":/a10f5:"y"::."Voronoi Diagram":WeSn -X5 -Y5 -P -K > $fout
awk '{print $1,$2}' $fS | psxy -R -J -W1 -Gred -Sc.1 -O -K >> $fout
awk '{print $1,$2}' $fE | psxy -R -J -A -m -Ccolor.cpt -W2 -O -K >> $fout
pwd | psxy -R -J -O >> $fout
