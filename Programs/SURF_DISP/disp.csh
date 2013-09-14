cd /home/tianye/for_weisen/Thick_Vs_global_merged
foreach file (`ls Thick_Vs_*.txt`)
set w_flag=`awk 'NR==1{if($2==0){print "-c 1"}else{print ""}}' $file`
awk '{if($2==0){print $1,1.5,$2,$2*1.73*0.32 + 0.77}else{print $1,$2*1.73,$2,$2*1.73*0.32 + 0.77}}' $file > /home/tianye/code/Programs/SURF_DISP/temp.txt
cd /home/tianye/code/Programs/SURF_DISP
SURF_DISP temp.txt $file R 1 1 8 15 1 $w_flag
mv $file'.R.grv' /home/tianye/for_weisen/Thick_Vs_global_merged_grv
rm -f $file'.R.phv' $file'.R'
cd /home/tianye/for_weisen/Thick_Vs_global_merged
end
