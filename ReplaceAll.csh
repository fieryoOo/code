set tdir = Programs
#set tdir = bin
set phome = `echo ${HOME} | sed s/'\/'/'\\\/'/g`

set flog = replace.log
rm -f $flog
#foreach file ( `find $tdir -type f \( -name "*\.c*" -o -name "*\.cpp*" -o -name "*\.C*" -o -name "*\.f*" -o -name "*\.csh*" \)` )
foreach file ( `find $tdir -type f` )
#foreach file (Script/FTAN/FTAN_SNR_ETHQ_fixed_bp.csh)
if( `grep '/home/tianye' $file | grep -v '/code' | grep -v '/Model' | grep -v '/MyLib' | grep -v '/home/tianye/for_' | grep -c -v '/Software' ` == 0 ) then
continue
endif
echo "Searching in "$file
set ftmp = $file'_replace_temp'
echo "#########################"$file"#########################" >> $flog
sed s/'\/home\/tianye\/Programs'/'\/home\/tianye\/code\/Programs'/g $file | sed s/'\/home\/tianye\/Script'/'\/home\/tianye\/code\/Script'/g | sed s/'\/home\/tianye\/sac'/'\/home\/tianye\/Software\/sac'/g | sed s/'\/home\/tianye\/code\/Programs\/pssac'/${phome}'\/usr\/bin'/g > $ftmp
#comm -23 $file $ftmp >> $flog
#echo "---------- to ----------" >> $flog
#comm -13 $file $ftmp >> $flog
sdiff -s -w 148 $file $ftmp >> $flog
echo "#\n#" >> $flog
#rm -f $ftmp
mv $ftmp $file
end
