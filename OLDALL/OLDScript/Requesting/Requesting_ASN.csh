#!/bin/csh
if ( $#argv != 5 )then
echo "usage: Requesting_ASN.csh [sta_list (as 'M12A TA')] [month_list (as '2008 6')] [label] [email] [name]"
exit
endif

#set sta_num=`more $argv[1] | wc -l`
#set lst_num=`echo $sta_num | awk '{print int(($1-1)/50)+1}'`
#set i=1
#while ( $i <= $lst_num )
#awk -v i=$i 'NR>(i-1)*50 && NR<=i*50' $argv[1] > 'temp_sta_lst_'$i
#@ i += 1
#end #while

foreach ymonth (`awk '{print $1"@"$2}' $argv[2]`)

set year=`echo $ymonth | cut -d@ -f1`
set month=`echo $ymonth | cut -d@ -f2`
 if( $month == 1 ) then
	set month_name =  JAN
	set day_num = 31
 endif
 if( $month == 2 ) then
	set month_name =  FEB
	set day_num = 28
        if( `echo $year | awk '{if(($1%4==0 && $1%100!=0) || $1%400==0){print 1}}'` ) set day_num = 29
 endif
 if( $month == 3 ) then
	set month_name =  MAR
	set day_num = 31
 endif
 if( $month == 4 ) then
	set month_name =  APR
	set day_num = 30
 endif
 if( $month == 5 ) then
	set month_name =  MAY
	set day_num = 31
 endif
 if( $month == 6 ) then
	set month_name =  JUN
	set day_num = 30
 endif
 if( $month == 7 ) then
	set month_name =  JUL
	set day_num = 31
 endif
 if( $month == 8 ) then
	set month_name =  AUG
	set day_num = 31
 endif
 if( $month == 9 ) then
	set month_name =  SEP
	set day_num = 30
 endif
 if( $month == 10 ) then
	set month_name =  OCT
	set day_num = 31
 endif
 if( $month == 11 ) then
	set month_name =  NOV
	set day_num = 30
 endif
 if( $month == 12 ) then
	set month_name =  DEC
	set day_num = 31
 endif

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


set day=1
while ( $day <= $day_num )
#set day=12
#while ( $day <= 12 )
echo ".NAME "$argv[5]"\n.INST CU\n.MAIL University of Colorado at Boulder\n.EMAIL "$argv[4]"\n.PHONE\n.FAX\n.MEDIA: Electronic (FTP)\n.ALTERNATE MEDIA: Electronic (FTP)\n.LABEL "$argv[3]"_"$year"."$month_name"."$day"\n.END" > requesting_email
if ( $day < $day_num ) then
set day2=`echo $day | awk '{print $1+1}'`
set temp=$year' '$month' '$day' 0 0 0 '$year' '$month' '$day2' 0 0 0 1 LHZ'
else
set temp=$year' '$month' '$day' 0 0 0 '$year2' '$month2' 1 0 0 0 1 LHZ'
endif
foreach sta_n (`awk '{print $1"@"$2}' $argv[1]`)
set sta=`echo $sta_n | cut -d@ -f1`
set net=`echo $sta_n | cut -d@ -f2`
echo $sta' '$net' '$temp >> requesting_email
end #sta_n
cat requesting_email | mail -s 'Requesting Data' breq_fast@iris.washington.edu
#set sta_b=`echo $i | awk '{print ($1-1)*50+1}'`
#set sta_e=`echo $i | awk '{print $1*50}'`
echo 'Request of '$year $month_name $day' sent...'
@ day += 1
end #while

end #ymonth

#rm -rf temp_sta_lst_*
