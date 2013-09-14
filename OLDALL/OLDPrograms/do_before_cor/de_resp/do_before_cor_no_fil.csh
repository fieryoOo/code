#foreach  year  (2007)
#set year = 2007
set year = 2008
set bp = 5to150_new
foreach month (AUG)
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

	cd $year$month
#   cp ~/Tianshan/Info/station.lst .	#cp ~/Tianshan/Info/station.lst .
#if ( 1 == 0 ) then
	rm event.dat
#a detail of station.lst needs to be modified. More specificlly, province name should be added before station name, like "ProvStaName", eg. HLBAQ.#### 
	\rm event.dat
#	foreach day1 ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  )
	foreach day1 ( 12 )
	echo "remove ft*\n"
	rm ${year}_${nmonth}_${day1}_0_0_0/ft_*	
	echo $nmonth $day1 > temp.dat
	
	awk -v ytp=$year '{  printf "%-4s %2s %2s 000000.00\n",ytp,$1, $2}' temp.dat  >> event.dat
	end

	\rm sac_db.out
#	/mtera/weisen/for_yong/CODE/CHINA_DATA/set_sacdb_China_v2 LHZ
#	/mtera/weisen/for_yong/CODE/CHINA_DATA/CUT_TRANS_Zcomp/cut_trans_RESP_CH 1000 83000
	/Users/jiayixie/progs/prog_CV/set_sacdb_China_v2 BHZ
	/Users/jiayixie/progs/prog_CV/CUT_TRANS_Zcomp/cut_trans_RESP_CH 1000 83000	
#	reverse stnm_reverse.txt
	
	rm -r $bp
#endif
if ( 1 == 0) then
	mkdir $bp
	cd $bp



	cp ../*.*  .

	foreach day ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  )
		set dirday =  ${year}_${nmonth}_${day}_0_0_0
		mkdir  $dirday
		cd $dirday
		cp ../../${dirday}/ft_* .
		ls ft_*SAC > temp.ft
		awk '{print "200 150 5 4 1 1  ", $1}' temp.ft > param.dat
#		/mtera/weisen/for_yong/CODE/FTA/filter4_f/filter4  param.dat
		/Users/jiayixie/progs/prog_CV/FTA/filter4_f/filter4  param.dat

#		/mtera/weisen/for_yong/CODE/FTA/white_outphamp/whiten_rej_phamp  param.dat
		/Users/jiayixie/progs/prog_CV/FTA/white_outphamp/whiten_rej_phamp  param.dat
		\rm ft_*.SAC
		\rm ft_*.SAC_bit
		cd ../
	end

	\rm sac_db.out
#	/home/zheng/progs/NOISE_CODA/SAC_FROM_SEED/set_sacdb.am
#	/mtera/weisen/for_yong/CODE/SAC_FROM_SEED/set_sacdb.am
##	/Users/jiayixie/progs/prog_CV/SAC_FROM_SEED/set_sacdb.am
##	\rm COR -fr
##	mkdir COR
#	/home/zheng/progs/NOISE_CODA/NEWCORR/justCOR  5000
#	/mtera/weisen/for_yong/CODE/NEWCORR/justCOR  5000
##	/Users/jiayixie/progs/prog_CV/NEWCORR/justCOR  5000
	cd ../../
end
endif
