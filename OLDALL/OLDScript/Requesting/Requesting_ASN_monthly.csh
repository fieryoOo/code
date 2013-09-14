#!/bin/csh
if ( $#argv != 5 )then
echo "usage: Requesting_ASN.csh [sta_list (as 'M12A TA')] [month_list (as '2008 6')] [label] [email] [name]"
exit
endif

set sta_num=`more $argv[1] | wc -l`
set lst_num=`echo $sta_num | awk '{print int(($1-1)/50)+1}'`
set i=1
while ( $i <= $lst_num )
awk -v i=$i 'NR>(i-1)*50 && NR<=i*50' $argv[1] > 'temp_sta_lst_'$i
@ i += 1
end #while

foreach ymonth (`awk '{print $1"@"$2}' $argv[2]`)

set year=`echo $ymonth | cut -d@ -f1`
set month=`echo $ymonth | cut -d@ -f2`
 if( $month == 1 )  set month_name =  JAN
 if( $month == 2 )  set month_name =  FEB
 if( $month == 3 )  set month_name =  MAR
 if( $month == 4 )  set month_name =  APR
 if( $month == 5 )  set month_name =  MAY
 if( $month == 6 )  set month_name =  JUN
 if( $month == 7 )  set month_name =  JUL
 if( $month == 8 )  set month_name =  AUG
 if( $month == 9 )  set month_name =  SEP
 if( $month == 10 ) set month_name =  OCT
 if( $month == 11 ) set month_name =  NOV
 if( $month == 12 ) set month_name =  DEC

if ( $month > 0 && $month < 12 )then
set year2=$year
set month2=`echo $month | awk '{print $1+1}'`
else if ( $month == 12 )then
set year2=`echo $year | awk '{print $1+1}'`
set month2=1
else
echo "Wrong month info!"
continue
endif

set i=1
while ( $i <= $lst_num )
echo ".NAME "$argv[5]"\n.INST CU\n.MAIL University of Colorado at Boulder\n.EMAIL "$argv[4]"\n.PHONE\n.FAX\n.MEDIA: Electronic (FTP)\n.ALTERNATE MEDIA: Electronic (FTP)\n.LABEL "$argv[3]"_"$year"."$month_name"\n.END" > requesting_email
foreach sta_n (`awk '{print $1"@"$2}' 'temp_sta_lst_'$i`)
set sta=`echo $sta_n | cut -d@ -f1`
set net=`echo $sta_n | cut -d@ -f2`
echo $sta' '$net' '$year' '$month' 1 0 0 0 '$year2' '$month2' 1 0 0 0 1 LHZ' >> requesting_email
end #sta_n
cat requesting_email | mail -s 'Requesting Data' breq_fast@iris.washington.edu
set sta_b=`echo $i | awk '{print ($1-1)*50+1}'`
set sta_e=`echo $i | awk '{print $1*50}'`
echo 'Request of '$year $month_name' for the '$sta_b' to '$sta_e'th station sent...'
@ i += 1
end #while

end #ymonth

rm -rf temp_sta_lst_*
