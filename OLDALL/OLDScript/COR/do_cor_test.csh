#!/bin/csh
set mainpath = /utera/tianye/data_test_method
set station_lst = ${mainpath}/station.lst

cd ${mainpath}
foreach year (2008)
#set year = 2008
#set bp = 5to150
foreach bp ( 5to150_onebit_whiten )
#foreach month ( JAN FEB MAR APR MAY JUN)
foreach month ( APR )
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
	/home/tianye/code/Programs/CORR_REC/set_sacdb.am ${station_lst} ${event_dat}
	\rm COR -fr
	mkdir COR
	#/home/tianye/code/Programs/CORR_REC/justCOR 3000
/home/tianye/code/Script/COR/NEWCORR/justCOR_part_station 3000
	cd ../../
end
end
end #bp
