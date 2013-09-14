gmtset HEADER_OFFSET -0.1
gmtset LABEL_OFFSET -0.1
gmtset HEADER_FONT_SIZE 14
gmtset LABEL_FONT_SIZE 10
gmtset ANNOT_FONT_SIZE_PRIMARY 8

#set Qm = `pwd | awk -F/ '{print $NF}'`
set Qm = QF
set dir0 = UNC
set dir1 = separate_${Qm}1s
set dir2 = separate_${Qm}Melt1s
set REG = -R3.9/4.6/10/50
set SCA = -JX1.5i/-2.5i
set Xs = 5
set clr = {30/30/30,blue,red}

set outf = Vs_comparison_${Qm}.ps
rm -f $outf
set N=0
foreach age (0 0.5 1 1.5 2 2.5 3 3.5)
#foreach age (0)
set Y = 0
if( $N == 4) then
set Shift = '-X-15 -Y-10'
else
set Shift = -X$Xs
endif
if( $age == 0 ) then
set f0 = Vs_obs/MC.1.0.0.mod
set f1 = ${dir1}/age0.txt
set f2 = ${dir2}/age0.txt
psbasemap $REG $SCA -Ba0.3f0.1:"Vs(km/sec)":/a10f5:"Depth (km)"::."$age Ma":WSne -X5 -Y12 -K >> $outf
awk '{print $2,$1,$3}' $f0 | psxy $REG $SCA -W4,$clr[1] -A -O -K >> $outf #-Ex.1c
else
set agetmp = `echo $age | awk '{printf "%.1f",$1}'`
set f0 = ${dir0}/${agetmp}.mod.out
set f1 = ${dir1}/age${age}.txt_sm
set f2 = ${dir2}/age${age}.txt_sm
psbasemap $REG $SCA -Ba0.3f0.1:"Vs(km/sec)":/a10f5:"Depth (km)"::."$age Ma":WSne $Shift -O -K >> $outf
awk '{print $1,$2+2*$3}' $f0 > temp1
awk '{print $1,$2-2*$3}' $f0 > temp2
~/code/Programs/Smoothing/Fast_Smoothing_curve temp1 0.5
~/code/Programs/Smoothing/Fast_Smoothing_curve temp2 0.5
awk 'BEGIN{i=0}{x[i++]=$0}END{for(i=i-1;i>=0;i--)print x[i]}' temp2_sm >> temp1_sm
awk '{print $2,$1}' temp1_sm | psxy $REG $SCA -W2,$clr[1] -Ggray -A -L -O -K >> $outf
endif
awk '{print $2,$1}' $f2 | psxy $REG $SCA -W4,$clr[3],- -A -O -K >> $outf
awk '{print $2,$1}' $f1 | psxy $REG $SCA -W4,$clr[2],- -A -O -K >> $outf
@ N ++
end #foreach age
rm -f temp1 temp2 temp1_sm temp2_sm
pwd | psxy $REG $SCA -O >> $outf
echo $outf
