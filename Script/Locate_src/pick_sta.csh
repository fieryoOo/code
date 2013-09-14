rm -f station_TA_exist.lst
rm -f station.lst

awk -F/ '{if($1!=sta){print $1; sta=$1}}' file_daynum.lst > temp.lst
foreach sta (`tail file_daynum.lst | awk -F_ '{print $3}' | awk -F. '{print $1}'`)
if( ! `awk -v sta=$sta '{if($1==sta){print 1; exit}}' temp.lst` ) echo $sta >> temp.lst
end

foreach sta (`more temp.lst`)
set location=`awk -v sta=$sta '{if($1==sta){print $2"@"$3; exit}}' station_TA.lst`
set long=`echo $location | cut -d@ -f1`
set lati=`echo $location | cut -d@ -f2`
echo $sta $long $lati >> station_TA_exist.lst
set dist1=`/home/tianye/code/Programs/DIST/get_dist 40 236 $lati $long`
set dist2=`/home/tianye/code/Programs/DIST/get_dist 50 260 $lati $long`
if (`echo $dist1 $dist2 | awk '{if( ($1<1050 && $1>200) || ($2<900 && $2>200)){print 1}else{print 0}}'`) echo $sta $long $lati >> station.lst
end

