#!/bin/csh
if($#argv != 3)then
echo "Usage: src_trvtime_input_path.csh [dir (with event.dat)] [sta_path.lst(put into dir/Predict_trvt)] [env_dir]"
exit
endif

set dir=$argv[1]
cd $dir
set sta_path=`echo $argv[2] | cut -d. -f1`
set per=12
set sta1=`awk 'NR==1 {print $1}' Predict_trvt/$argv[2]`
set loc1=`awk 'NR==1 {print $3,$2}' Predict_trvt/$argv[2]`
###event T2457###
#set long_e=236
#set lati_e=40
###event T2870###
#set long_e=244
#set lati_e=45

rm -rf Predict_trvt/$sta_path
mkdir Predict_trvt/$sta_path
rm -f 'Predict_trvt/'$sta_path'/dist_dt_temp'
set flag=0
foreach location ( 236@40 239@45 290@55 ) #T2457 T2865 T3726
#foreach location ( 236@40 244@45 244@47 290@55 ) #T2457 T2870 T3032 T3726
  set long_e=`echo $location | cut -d@ -f1`
  set lati_e=`echo $location | cut -d@ -f2`

  set event=`awk -v long=$long_e -v lati=$lati_e '$2==long && $3==lati {print $1}' event.dat`
echo $event
  set loce=`echo $lati_e $long_e`
  set dist1=`/home/tianye/code/Programs/DIST/get_dist $loce $loc1`
  set v1=`awk -v per=$per '$1==per {print $2}' '/home/tianye/Model/Cst2_Spr/vel/Pre_T/'$event'/'$event'_'$sta1'.dat'`
  set t1=`echo $dist1 $v1 | awk '{print $1/$2}'`

#  rm -rf Predict_trvt/$sta_path
#  mkdir Predict_trvt/$sta_path
#  rm -f 'Predict_trvt/'$sta_path'/dist_dt_temp'
  echo 'src_dt_'$event'_'$sta1'.ps'
  foreach info2 (`awk 'NR>1 {print $1"@"$2"@"$3}' Predict_trvt/$argv[2]`)
    set sta2=`echo $info2 | cut -d@ -f1`
    if( $flag ) goto next
    set file=$argv[3]'/'$sta1'/COR_'$sta1'_'$sta2'.SAC.env'
    if( -e $file ) then
      cp $file Predict_trvt/$sta_path
    else
      set file=$argv[3]'/'$sta2'/COR_'$sta2'_'$sta1'.SAC.env'
      if( -e $file ) then
        set out_env='Predict_trvt/'$sta_path'/COR_'$sta1'_'$sta2'.SAC.env'
echo $out_env
sac << END
r $file
reverse
w $out_env
quit
END
      endif
    endif
next:
    if( ! -e 'Predict_trvt/'$sta_path'/COR_'$sta1'_'$sta2'.SAC.env')then
      echo "no env found between "$sta1" and "$sta2
      continue
    endif
    set loc2=`echo $info2 | awk -F@ '{print $3,$2}'`
    set dist2=`/home/tianye/code/Programs/DIST/get_dist $loce $loc2`
    set v2=`awk -v per=$per '$1==per {print $2}' '/home/tianye/Model/Cst2_Spr/vel/Pre_T/'$event'/'$event'_'$sta2'.dat'`
    set t2=`echo $dist2 $v2 | awk '{print $1/$2}'`
    set dt=`echo $t1 $t2 | awk '{print $2-$1}'`
    set sta_dist=`/home/tianye/code/Programs/DIST/get_dist $loc1 $loc2`
    echo $sta_dist $dt 'COR_'$sta1'_'$sta2'.SAC.env' >> 'Predict_trvt/'$sta_path'/dist_dt_temp'
  end #foreach info2
  set flag=1

  cd 'Predict_trvt/'$sta_path
  sort -gr dist_dt_temp > dist_dt_temp2
  mv dist_dt_temp2 dist_dt_temp
  set dmax=`head -1 dist_dt_temp | awk '{print $1+100}'`
  set JXY = -JX10/15
  set REG = -R-1000/1000/0/$dmax

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
  cd ../..
end #location
