set Npath=`more file_daynum.lst | wc -l`
rm -f n.env.lst
set sta='AAAA'
foreach info (`awk '{print $1}' file_daynum.lst`)
#set dnum=`echo $info | cut -d@ -f2`
#if ( $dnum < 50 )continue
set sta1=`echo $info | cut -d/ -f1`
set sta2=`echo $info | cut -d_ -f3 | cut -d. -f1`
if ( $sta1 == $sta ) then
  if ( ! $flag ) if(! `awk -v sta2=$sta2 '{if($1==sta2){print 1;exit}}' station.lst`) continue
else
  set sta=$sta1
  set flag=`awk -v sta1=$sta1 '{if($1==sta1){print 1;exit}}' station.lst`
  if( ! $flag )if(! `awk -v sta2=$sta2 '{if($1==sta2){print 1;exit}}' station.lst`) continue
endif
set file=`echo $info | cut -d@ -f1`
#set N=`echo $info | cut -d@ -f2`
echo $file" 12 "$sta1 $sta2 >> n.env.lst
#echo $N"/"$Npath finished...
end
