set dir='/home/tianye/Model/GR_PH_maps/Vs_Cst2_Shapiro_merged'
cd $dir
foreach file (`ls Merged_Vs_*.txt`)
#foreach file ( Merged_Vs_102_-56.txt )
#foreach file ( Merged_Vs_114_-62.txt )
set w_flag=`awk 'NR==1{if($2==0){print "-c 1"}else{print ""}}' $file`
awk '{if($2==0){print $1,1.5,$2,$2*1.73*0.32 + 0.77}else{print $1,$2*1.73,$2,$2*1.73*0.32 + 0.77}}' $file > /home/tianye/code/Programs/SURF_DISP.dir/temp.txt
cd /home/tianye/code/Programs/SURF_DISP.dir

SURF_DISP temp.txt $file R 1 1 5 19 1 $w_flag
mv $file'.R.grv' gtemp1
mv $file'.R.phv' ptemp1
SURF_DISP temp.txt $file R 1 1 5 100 5 $w_flag
mv $file'.R.grv' gtemp2
mv $file'.R.phv' ptemp2

awk '$1>=20' gtemp2 >> gtemp1
awk '$1>=20' ptemp2 >> ptemp1
awk 'NF>0' gtemp1 > $dir'/'$file'.R.grv'
awk 'NF>0' ptemp1 > $dir'/'$file'.R.phv'

rm -f gtemp? ptemp? $file'.R'
cd $dir
end
