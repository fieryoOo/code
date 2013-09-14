#!/bin/csh
set mainpath = /utera/tianye/US_ASN_amp
set station_lst = ${mainpath}/station.lst

cd ${mainpath}
foreach year (2009)
#set year = 2008
set bp = 5to150
#foreach month ( JAN FEB MAR APR MAY JUN)
foreach month ( JAN )
	if( $month == JAN ) @ nmonth =  1
	if( $month == FEB ) @ nmonth =  2
	if( $month == MAR ) @ nmonth =  3
	if( $month == APR ) @ nmonth =  4
	if( $month == MAY ) @ nmonth =  5
	if( $month == JUN ) @ nmonth =  6
	if( $month == JUL ) @ nmonth =  7
	if( $month == AUG ) @ nmonth =  8
	if( $month == SEP ) @ nmonth =  9
	if( $month == OCT ) @ nmonth = 10
	if( $month == NOV ) @ nmonth = 11
	if( $month == DEC ) @ nmonth = 12

	cd $year'.'$month

	cd $bp
	set event_dat = ${mainpath}/$year'.'$month/event.dat

    \rm sac_db.out
echo aaa
	/home/tianye/code/Programs/CORR_REC/set_sacdb.am ${station_lst} ${event_dat}
	\rm COR -fr
	mkdir COR
echo AAA
	/home/tianye/code/Programs/CORR_REC/justCOR 3000
echo BBB
	cd ../../
end
end
