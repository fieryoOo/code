#!/bin/csh
set bp = 5to150
#set year = 2007
set label1 = OBS_7D
set label2 = OBS_US

set mainpath = /media/WORK/tianye/ASN_OBS
set station_lst = ${mainpath}/station.lst
cd $mainpath
foreach year ( 2011 2012 )
foreach month ( JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC )
#foreach month ( NOV )

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

mkdir -p $year'.'$month
set seed_flag=0
#foreach label ( $label1 $label2 )
foreach label ( $label1 )
rm -f $year'.'$month'/input_ev_seed'
###preparing SEED.lst
foreach day1 ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 )
#foreach day1 ( 3 )
echo 'SEED/'$label'_'$year'.'$month'.'$day1
mv SEED/$label'_'$year'.'$month'.'$day1'.'*'.seed'  SEED/$label'_'$year'.'$month'.'$day1
pwd
echo SEED/$label'_'$year'.'$month'.'$day1
if ( -e SEED/$label'_'$year'.'$month'.'$day1 ) then
echo " PDE   ${year}    $nmonth $day1  0000000000000     63.52 -147.44  11 8.50 Ms GS   9C.G F.S."  >> $year'.'$month'/input_ev_seed'
echo ../SEED/$label'_'$year'.'$month'.'$day1 >> $year'.'$month'/input_ev_seed'
endif
end #foreach day1

if( ! -e $year'.'$month'/input_ev_seed' ) continue
set seed_flag = 1

cd $year'.'$month
cp $station_lst ./station.lst
/home/tianye/code/Programs/do_before_cor/de_resp/sac_from_seed_RESP LHZ
mv event_station.tbl event_station_${label}.tbl
#/home/tianye/code/Programs/do_before_cor/de_resp/cut_trans_RESP 1000 84000
cd ..
end #foreach label
if( ! $seed_flag ) then
echo "No data for "$month"/"$year". Skipped."
rm -rf $year'.'$month
continue
endif
end
end
exit
###Test data existence
if( ! $seed_flag ) then
echo "No data for "$month"/"$year". Skipped."
rm -rf $year'.'$month
continue
endif
cd $year'.'$month

echo "Removing the original SAC and RESP files..."
foreach day ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 )
#foreach day (3)
rm -f ${year}.${month}.${day}/${year}.${month}.${day}.*.LHZ.SAC
rm -f ${year}.${month}.${day}/RESP.*.LHZ
end

rm -rf $bp
mkdir $bp
cd $bp

foreach day ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  )
#foreach day ( 3 )
set dirday =  ${year}.${month}.${day}
mkdir -p $dirday
cd $dirday
mv ../../${dirday}/ft_*.LHZ.SAC .
if( `ls ../../${dirday} | wc -l` == 0 )rm -rf ../../${dirday}
ls ft_* > temp.ft
awk '{print "60 50 10 7.5 1 1  ", $1}' temp.ft > param.dat
/home/tianye/code/Programs/do_before_cor/filter4/filter4 param.dat 1
awk '{print "100 80 4 3 1 1  ", $1, $1"_10.0_50.0"}' temp.ft > param.dat
/home/tianye/code/Programs/do_before_cor/whiten_outphamp/whiten_rej_phamp_no_fil  param.dat
rm -f param.dat temp.ft ft_*.SAC ft_*.SAC_10.0_50.0
cd ../
set nfile=`ls $dirday | wc -l`
if( $nfile == 0 )rm -rf $dirday
end

	rm -f event.dat
        foreach day1 ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  )
        echo $nmonth $day1 > temp.dat
        awk -v ytp=$year '{  printf "%-4s %2s %2s 000000.00\n",ytp,$1, $2}' temp.dat  >> event.dat
        end #foreach day1

/home/tianye/code/Programs/CORR_REC/set_sacdb.am ${station_lst} event.dat
exit
mkdir -p COR

/home/tianye/code/Programs/CORR_REC/justCOR 3000 0

cd COR
foreach dir (`ls -d */`)
if( `ls $dir | wc -l` == 0 )rm -rf $dir
end
cd ..

rm -f temp.dat sac_db.out event.dat
cd ..
rm -f from_seed input_ev_seed list_sac1 sac_bp_respcor sac_db.out
cd ..

end #foreach month

end #foreach year
