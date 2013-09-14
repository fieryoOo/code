set dir='/home/tianye/data_locate_src'
set per=12
#set sta1='V11A'
foreach sta1 ( A05A E06A V11A )
###event T2457###
#set long_e=236
#set lati_e=40
###event T2870###
#set long_e=244
#set lati_e=45

foreach location ( 236@40 244@45 ) #T2457 T2870
set long_e=`echo $location | cut -d@ -f1`
set lati_e=`echo $location | cut -d@ -f2`

set JXY = -JX10/15
set REG = -R-1000/1000/0/3000

cd $dir

set loc1=`awk -v sta=$sta1 '$1==sta' station.lst | awk '{print $3,$2}'`
set event=`awk -v long=$long_e -v lati=$lati_e '$2==long && $3==lati {print $1}' event.dat`
set loce=`echo $lati_e $long_e`
set dist1=`/home/tianye/code/Programs/DIST/get_dist $loce $loc1`
set v1=`awk -v per=$per '$1==per {print $2}' '/home/tianye/Model/Cst2_Spr/vel/Pre_T/'$event'/'$event'_'$sta1'.dat'`
set t1=`echo $dist1 $v1 | awk '{print $1/$2}'`

foreach file (`ls *'/COR_'*'_'$sta1'.SAC.env'`)
set sta2=`echo $file | cut -d/ -f1`
set out_env=$sta1'/COR_'$sta1'_'$sta2'.SAC.env'
if( -e $out_env ) continue
sac << END
r $file
reverse
w $out_env
quit
END
end

rm -f dist_dt_temp
echo 'src_dt_'$event'_'$sta1'.ps'

foreach file (`ls $sta1'/COR_'$sta1*'.SAC.env'`)
set sta2=`echo $file | cut -d_ -f3 | cut -d. -f1`
set loc2=`awk -v sta=$sta2 '$1==sta' station.lst | awk '{print $3,$2}'`
set dist2=`/home/tianye/code/Programs/DIST/get_dist $loce $loc2`
set v2=`awk -v per=$per '$1==per {print $2}' '/home/tianye/Model/Cst2_Spr/vel/Pre_T/'$event'/'$event'_'$sta2'.dat'`
set t2=`echo $dist2 $v2 | awk '{print $1/$2}'`
set dt=`echo $t1 $t2 | awk '{print $2-$1}'`
set sta_dist=`/home/tianye/code/Programs/DIST/get_dist $loc1 $loc2`
echo $sta_dist $dt $file >> dist_dt_temp
end

sort -gr dist_dt_temp > dist_dt_temp2
mv dist_dt_temp2 dist_dt_temp

if ( ! -e src_dt_'$sta1'_origin2.ps )then
awk '{print $3}' dist_dt_temp | pssac $REG ${JXY} -Ba500f100/a500f100WeSn:.$event'_'$sta1':' -M0.3 -Ekt-1 -Q > 'src_dt_'$sta1'_origin.ps'
awk '{print $3}' dist_dt_temp | pssac $REG ${JXY} -Ba500f100/a500f100WeSn:.$event'_'$sta1':' -M0.3 -Ekt-1 -Q -G155/155/255/0.00006 > 'src_dt_'$sta1'_origin2.ps'
endif

awk '{print $3}' dist_dt_temp | pssac ${REG} ${JXY} -Ba500f100/a500f100WeSn:.$event'_'$sta1':' -M0.3 -Q -Ekt-1 -G155/155/255/0.00006 -K > 'src_dt_'$event'_'$sta1'.ps'

foreach dot (`awk '{print $2"@"$1+100}' dist_dt_temp`)
set dt=`echo $dot | cut -d@ -f1`
set dist=`echo $dot | cut -d@ -f2`
echo $dt $dist | psxy $REG ${JXY} -Sc.1 -W.5 -G255/0/0 -K -O >> 'src_dt_'$event'_'$sta1'.ps'
end

pwd | psxy $REG -J -O >> 'src_dt_'$event'_'$sta1'.ps'

rm -f dist_dt_temp

end #location
end #sta1
