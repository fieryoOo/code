set i=0
set lati=10
while ( $lati <= 70 )
echo $lati

set long=210
while ( $long <= 290 )
#set flag=`echo $lati $long | awk '{if(($1>-1.16129*$2+289.871)&&(($1<-1.66667*$2+470)||($1<45))){print 1}else{print 0}}'`
#if ( $flag == 0 ) then
#@ long += 1
#continue
#endif
@ i += 1
echo $i'\tT'$i'\t'$long'\t'$lati >> event.dat
@ long += 1
end

@ lati += 1
end
