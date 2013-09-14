#awk '{print $1}' file_daynum_30.lst | awk -F_ '{print $2,$3,$0}' | awk -F. '{print $1,$2".SAC"}' | awk '{print $4".env",12,$1,$2}' > n.env.lst
/home/tianye/code/Programs/Locate_Source/CODE/get_source_YONG n.env.lst station.lst event.dat /data/ulisse/tianye/Model/Crust2_Shapiro/vel_map/Prediction_T

rm -rf result_T
mkdir result_T
mv 12_*.result result_T
cd result_T
ls 12_*.result > result.lst
set f0 = `head -n1 result.lst`
#python /home/tianye/code/Programs/Locate_Source/CODE/do_normalization.py $f0 > n_$f0

cp $f0 mtmp.txt
#foreach file ( `awk '{if (NR<=500 && NR >=1) print $1}' result.lst` )
foreach file ( `awk '{if (NR >=2) print $1}' result.lst` )
echo $file
#python /home/tianye/code/Programs/Locate_Source/CODE/do_normalization.py $file > n_$file
python /home/tianye/code/Programs/Locate_Source/Precursor_Weisen/cvadd.py $file mtmp.txt  > mtmp
mv mtmp mtmp.txt

end

