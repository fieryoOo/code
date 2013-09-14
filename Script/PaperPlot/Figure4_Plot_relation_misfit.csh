### Grid search result
gmtset HEADER_OFFSET 0
gmtset LABEL_OFFSET 0
#gmtset LABEL_FONT_SIZE 18
#gmtset HEADER_FONT_SIZE 20
gmtset ANNOT_FONT_SIZE_PRIMARY 14
gmtset LABEL_FONT 4
gmtset ANNOT_FONT 4
gmtset HEADER_FONT 4
set outf = 'Figure2_Relation_Misfit.ps'
set Xs = 8.5
set clr = {red,orange,forestgreen,blue,black}
####################
### Plot Group Relations ###
set SCA = -JX3i/2.5i
set REG = -R/0.5/4.0/0.5/4.0
psbasemap $REG $SCA -Ba1f0.5:"Age (Ma)":/a0.5f0.1:"Velocity (km/sec)":/a0.5f0.1:."Group-Age Relationships":WSne -X3 -Y8 -K > $outf
set i = 1
foreach per (7.0 8.0 10.0 15.0)
set x = 0
rm -f plot.tmp
while ( `echo $x | awk '{if($1<=5)print 1; else print 0}'` )
awk '$1=='${per}'{print '${x}',$9+$10*('${x}'**0.5)+$11*'${x}'}' misfit/relation_all_g_avg_per.txt >> plot.tmp
set x = `echo $x | awk '{print $1+0.1}'`
end #while x
#psxy plot.tmp $REG $SCA -S-0.05 -W4,$clr[$i] -A -O -K >> $outf
psxy plot.tmp $REG $SCA -W4,$clr[$i],- -A -O -K >> $outf
awk '{print $3,$2}' Age_Grv/Age_Grv_$per'sec_GS_1' | psxy $REG $SCA -Sc0.1 -W1,black -G$clr[$i] -A -O -K >> $outf
#awk '{print $1,$2}' Age_Grv/Age_Grv_$per'sec_wavg_1' | psxy $REG $SCA -Sc0.1 -W1,black -G$clr[$i] -A -O -K >> $outf
@ i ++
end #per
pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $outf
-1.5 12.5 20 0.0 4 LT (a)
EOF
####
####################
### Plot Phase Relations ### shift 1
set REG = -R/0.5/4.0/2/4.0
psbasemap $REG $SCA -Ba1f0.5:"Age (Ma)":/a0.5f0.1:"Velocity (km/sec)":/a0.5f0.1:."Phase-Age Relationships":WSne -X10.5 -K -O >> $outf
set i = 1
foreach per (7.0 8.0 10.0 15.0)
set x = 0
rm -f plot.tmp
while ( `echo $x | awk '{if($1<=5)print 1; else print 0}'` )
awk '$1=='${per}'{print '${x}',$9+$10*('${x}'**0.5)+$11*'${x}'}' misfit/relation_all_p_avg_per.txt >> plot.tmp
set x = `echo $x | awk '{print $1+0.1}'`
end #while x
psxy plot.tmp $REG $SCA -W4,$clr[$i] -A -O -K >> $outf
awk '{print $3,$2}' Age_Phv/Age_Phv_$per'sec_GS_1' | psxy $REG $SCA -Sc0.1 -W1,black -G$clr[$i] -A -O -K >> $outf
#awk '{print $1,$2}' Age_Phv/Age_Phv_$per'sec_wavg_1' | psxy $REG $SCA -Sc0.1 -W1,black -G$clr[$i] -A -O -K >> $outf
@ i ++
end #per
pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $outf
-1.5 12.5 20 0.0 4 LT (b)
EOF
####
### add legend ###
pslegend -R0/10/0/5 -D11/2.7/3/3/BL $SCA -O -K -G255 << EOF >> $outf
L 15 4 l group vel:
S .2i - .3i - 2p,$clr[4],4:0p .5i  15 sec
S .2i - .3i - 2p,$clr[3],4:0p .5i  10 sec
S .2i - .3i - 2p,$clr[2],4:0p .6i  8 sec
S .2i - .3i - 2p,$clr[1],4:0p .6i  7 sec
EOF
pslegend -R0/10/0/5 -D11/0.4/3/3/BL $SCA -O -K -G255 << EOF >> $outf
L 15 4 l phase vel:
S .2i - .3i - 2p,$clr[4] .5i  15 sec
S .2i - .3i - 2p,$clr[3] .5i  10 sec
S .2i - .3i - 2p,$clr[2] .6i   8 sec
S .2i - .3i - 2p,$clr[1] .6i   7 sec
EOF
#G 1.5i
#N 1
#D 0 1p
pwd | psxy -R -J -O >> $outf
echo $outf
exit
#######################
set REG = -R-30/30/0/50
set SCA = -JX2.5i/2.5i
set age = 0.5
### Plot 1st per of age-ref misfit ###
set per = 7.0
   ### misfit histogram # shift 3
#set REG = -R-15/15/0/50
set wd = 2
set meang1 = `awk '$1=="g" && $2=='$per' && $3=='$age'{printf "%.3g",$5}' misfit/mean_std_ref.txt`
set stdg1 = `awk '$1=="g" && $2=='$per' && $3=='$age'{printf "%.3g",$7}' misfit/mean_std_ref.txt`
set meanp1 = `awk '$1=="p" && $2=='$per' && $3=='$age'{printf "%.3g",$5}' misfit/mean_std_ref.txt`
set stdp1 = `awk '$1=="p" && $2=='$per' && $3=='$age'{printf "%.3g",$7}' misfit/mean_std_ref.txt`
#psbasemap $REG $SCA -Ba5f1:"Vel Misfit (%)":/a10f5:"Percentage (%)":/a10f5:."$per sec":WSne -X$Xs -Y4 -O -K >> $outf
#awk '{print $2}' 'misfit/Phv_'$per'sec_T_misfit_reg_age_'${age}'.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,blue -Gblue -O -K >> $outf
#awk '{print $2}' 'misfit/Grv_'$per'sec_T_misfit_reg_age_'${age}'.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,red -O -K >> $outf
#pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $outf
#0.8 9.5 10 0.0 4 LT  group_mean: $meang1  group_std: $stdg1
#0.8 8.8 10 0.0 4 LT  phase_mean: $meanp1  phase_std: $stdp1
#EOF
####################
### Plot 2nd per of age-ref misfit ###
set per = 15.0
   ### misfit histogram
#set REG = -R-20/20/0/50
set wd = 1
set meang2 = `awk '$1=="g" && $2=='$per' && $3=='$age'{printf "%.3g",$5}' misfit/mean_std_ref.txt`
set stdg2 = `awk '$1=="g" && $2=='$per' && $3=='$age'{printf "%.3g",$7}' misfit/mean_std_ref.txt`
set meanp2 = `awk '$1=="p" && $2=='$per' && $3=='$age'{printf "%.3g",$5}' misfit/mean_std_ref.txt`
set stdp2 = `awk '$1=="p" && $2=='$per' && $3=='$age'{printf "%.3g",$7}' misfit/mean_std_ref.txt`
#psbasemap $REG $SCA -Ba5f1:"Vel Misfit (%)":/a10f5:"Percentage (%)":/a10f5:."$per sec":WSne -Y-9 -O -K >> $outf
#awk '{print $2}' 'misfit/Phv_'$per'sec_T_misfit_reg_age_'${age}'.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,blue -Gblue -O -K >> $outf
#awk '{print $2}' 'misfit/Grv_'$per'sec_T_misfit_reg_age_'${age}'.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,red -O -K >> $outf
#pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $outf
#0.8 9.5 10 0.0 4 LT  group_mean: $meang2  group_std: $stdg2
#0.8 8.8 10 0.0 4 LT  phase_mean: $meanp2  phase_std: $stdp2
#EOF
####################
### Plot 1st per ###
set REG = -R/-30/30/0/50
set per = 7.0
   ### misfit histogram # shift 2
set wd = 2
set meang = `awk '$1=="g" && $2=='$per'{printf "%.3f",$5}' misfit/mean_std.txt`
set stdg = `awk '$1=="g" && $2=='$per'{printf "%.3f",$7}' misfit/mean_std.txt`
set meanp = `awk '$1=="p" && $2=='$per'{printf "%.3f",$5}' misfit/mean_std.txt`
set stdp = `awk '$1=="p" && $2=='$per'{printf "%.3f",$7}' misfit/mean_std.txt`
set var_red_g = `echo $stdg $stdg1 | awk '{printf "%.0f",($2**2-$1**2)/$2**2*100}'`
set var_red_p = `echo $stdp $stdp1 | awk '{printf "%.0f",($2**2-$1**2)/$2**2*100}'`
psbasemap $REG $SCA -Ba10f2:"Phase Vel Misfit (%)":/a10f5:"Percentage (%)":/a10f5:."Misfit  $per sec":WSne -X-10 -Y-9.2 -O -K >> $outf
awk '{print $2}' 'misfit/Phv_'$per'sec_T_misfit.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,blue -Gblue -O -K >> $outf
#awk '{print $2}' 'misfit/Grv_'$per'sec_T_misfit.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,red -O -K >> $outf
awk '{print $2}' 'misfit/Phv_'$per'sec_T_misfit_reg_age_'${age}'.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,red -O -K >> $outf
pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $outf
0.6 9.7 10 0.0 4 LT         Mean  RMS                      Var.reduction
0.35 9.0 10 0.0 4 LT  $meanp1     $stdp1       (0.5 Ma)
0.6 8.3 10 0.0 4 LT  $meanp   $stdp    (age-dep)
7.8 8.8 10 0.0 4 LT  $var_red_p%
-1 11.8 20 0.0 4 LT (b)
EOF
#1.7 9.0 10 0.0 4 LT  $meang       $stdg         $var_red_g%    (group)
####################
### Plot 2nd per ###
set per = 15.0
   ### misfit histogram
set REG = -R-10/10/0/55
set wd = 1
set meang = `awk '$1=="g" && $2=='$per'{printf "%.3f",$5}' misfit/mean_std.txt`
set stdg = `awk '$1=="g" && $2=='$per'{printf "%.3f",$7}' misfit/mean_std.txt`
set meanp = `awk '$1=="p" && $2=='$per'{printf "%.3f",$5}' misfit/mean_std.txt`
set stdp = `awk '$1=="p" && $2=='$per'{printf "%.3f",$7}' misfit/mean_std.txt`
set var_red_g = `echo $stdg $stdg2 | awk '{printf "%.0f",($2**2-$1**2)/$2**2*100}'`
set var_red_p = `echo $stdp $stdp2 | awk '{printf "%.0f",($2**2-$1**2)/$2**2*100}'`
psbasemap $REG $SCA -Ba5f1:"Phase Vel Misfit (%)":/a10f5:"Percentage (%)":/a10f5:."Misfit  $per sec":WSne -X$Xs -O -K >> $outf
awk '{print $2}' 'misfit/Phv_'$per'sec_T_misfit.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,blue -Gblue -O -K >> $outf
#awk '{print $2}' 'misfit/Grv_'$per'sec_T_misfit.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,red -O -K >> $outf
awk '{print $2}' 'misfit/Phv_'$per'sec_T_misfit_reg_age_'${age}'.txt' | pshistogram $REG $SCA -Z1 -W$wd -L1.0p,red -O -K >> $outf
pstext -R0/10/0/10 $SCA -V -O -K -N << EOF >>  $outf
0.6 9.7 10 0.0 4 LT         Mean  RMS                      Var.reduction
0.6 9.0 10 0.0 4 LT  $meanp2     $stdp2       (0.5 Ma)
0.6 8.3 10 0.0 4 LT  $meanp   $stdp    (age-dep)
7.8 8.8 10 0.0 4 LT  $var_red_p%
-1 11.8 20 0.0 4 LT (c)
EOF
#1.7 9.0 10 0.0 4 LT  $meang       $stdg         $var_red_g%    (group)
### ending ###
pwd | psxy $REG $SCA -O >> $outf
echo $outf
