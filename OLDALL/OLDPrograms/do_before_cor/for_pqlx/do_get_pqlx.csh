#foreach  year  (2007)
#set year = 2007
set year = 2007
set fileroot = /Users/jiayixie/Tianshan/Data_merged
set DATAROOT = /Users/jiayixie/Tianshan/Pqlx
set mainpass = `pwd`
set netwk = XJ
mkdir ${DATAROOT}/${netwk}

foreach month (OCT)
#foreach month (  SEP OCT NOV DEC )
#foreach month ( JAN FEB MAR APR MAY JUN JUL AUG SEP )
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

	cd ${fileroot}/$year$month
	cp ~/Tianshan/Info/station.lst .	#cp ~/Tianshan/Info/station.lst .
#if ( 1 == 0 ) then
	rm event.dat
#a detail of station.lst needs to be modified. More specificlly, province name should be added before station name, like "ProvStaName", eg. HLBAQ.#### 
#	foreach day1 (1 2)
	foreach day1 ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  )
	
	echo $nmonth $day1 > temp.dat
	
	awk -v ytp=$year '{  printf "%-4s %2s %2s 000000.00\n",ytp,$1, $2}' temp.dat  >> event.dat
	end#day1

	\rm sac_db_throw1.out
#	/Users/jiayixie/progs/jy/set_sacdb_China_throw1  BHZ /Users/jiayixie/Test/Get_SAC
	/Users/jiayixie/progs/jy/set_sacdb_China_throw1  BHZ ${DATAROOT}
	echo "=================begin to cut============\n"
#	/Users/jiayixie/progs/jy/cut_to_day_throw_1sig 1000 83000 /Users/jiayixie/Test/Get_SAC
	/Users/jiayixie/progs/jy/cut_to_day_throw_1sig 1000 83000 ${DATAROOT}
    cd ../
#endif

	end#month
#====================================
#pwd
#mkdir ${DATAROOT}/${netwk}/RESP
#cd ${DATAROOT}/${netwk}/RESP
#ls ~/Tianshan/RESP/RESP*BHZ | awk '{x=$0;split(x,a,"XJ");print "ln -s",$0   ,"'${DATAROOT}'/'${netwk}'/RESP/RESP.XJ." a[2]}' | sh

