rm -f TA.site
foreach info (`awk -F'\t' '{print $2,$3,$4,$5,$6,$7/1000.,$8}' station.info | awk '{print $1"@"$2"@"$4"@"$6"@"$7"@"$8"@"$9"@"$10"@"$11"@"$12"@"$13"@"$14"@"$15"@"$16"@"$17}'`)
echo "info: "$info
set be=`echo $info | cut -d@ -f1`
set t1=`echo $info | cut -d@ -f2`
set t2=`echo $info | cut -d@ -f3`
set mi=`echo $info | cut -d@ -f4-6 | sed s/'@'/' '/g`
set en=`echo $info | cut -d@ -f7-15 | sed s/'@'/' '/g`
set y=`echo $t1 | cut -d- -f1`
set m=`echo $t1 | cut -d- -f2`
set d=`echo $t1 | cut -d- -f3`
set jd=`~/code/Programs/mini_tools/calc_jday $y $m $d`
set tm1=`echo $y$jd`
set y=`echo $t2 | cut -d- -f1`
set m=`echo $t2 | cut -d- -f2`
set d=`echo $t2 | cut -d- -f3`
set jd=`~/code/Programs/mini_tools/calc_jday $y $m $d`
set tm2=`echo $y$jd`
echo $be $tm1 $tm2 $mi $en | awk '{printf "%-6s %8s %8s %9.4f %9.4f %9.4f@%s %s %s %s %s %s %s %s %s",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | awk -F@ '{printf "%s %-50s %-4s %-6s %9.4f %9.4f 09/13/12 17:37:44\n",$1,$2,"-","-",0,0}' >> TA.site
end

awk '{printf "%-6s %-8s %8s %8s %8s %-4s %9.4f %6.1f %6.1f %-50s 09/13/12 17:37:44\n%-6s %-8s %8s %8s %8s %-4s %9.4f %6.1f %6.1f %-50s 09/13/12 17:37:44\n%-6s %-8s %8s %8s %8s %-4s %9.4f %6.1f %6.1f %-50s 09/13/12 17:37:44\n",$1,"LHZ",$2,"-1",$3,"n",0,0,0,"-",$1,"LHE",$2,"-1",$3,"n",0,90,90,"-",$1,"LHN",$2,"-1",$3,"n",0,0,90,"-"}' TA.site > TA.sitechan
