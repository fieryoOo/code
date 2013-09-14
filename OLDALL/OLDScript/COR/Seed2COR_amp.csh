set bp = 5to150
#set year = 2007
set label = ASN_OBS2
set mainpath = /media/WORK/tianye/ASN_OBS
set station_lst = ${mainpath}/station.lst
cd $mainpath
foreach year ( 2011 2012 )
foreach month ( JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC )
#foreach month ( JAN )

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
foreach day1 ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 )
#foreach day1 ( 1 )
mv SEED/$label'_'$year'.'$month'.'$day1'.'*'.seed'  SEED/$label'_'$year'.'$month'.'$day1

if ( ! -e SEED/$label'_'$year'.'$month'.'$day1 ) continue
if ( $day1 == 1 ) then
echo " PDE   ${year}    $nmonth $day1  0000000000000     63.52 -147.44  11 8.50 Ms GS   9C.G F.S."  > $year'.'$month'/input_ev_seed'
else
echo " PDE   ${year}    $nmonth $day1  0000000000000     63.52 -147.44  11 8.50 Ms GS   9C.G F.S."  >> $year'.'$month'/input_ev_seed'
endif

echo ../SEED/$label'_'$year'.'$month'.'$day1 >> $year'.'$month'/input_ev_seed'

end

cd $year'.'$month
cp $station_lst ./station.lst

/home/tianye/code/Programs/do_before_cor/de_resp/sac_from_seed_RESP LHZ
/home/tianye/code/Programs/do_before_cor/de_resp/cut_trans_RESP 1000 84000
echo "Removing the original SAC and RESP files..."
foreach day ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  )
rm -f ${year}.${month}.${day}/${year}.${month}.${day}.*.LHZ.SAC
rm -f ${year}.${month}.${day}/RESP.*.LHZ
end

rm -rf $bp
mkdir $bp
cd $bp

foreach day ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  )
#foreach day ( 1 )
set dirday =  ${year}.${month}.${day}
mkdir -p $dirday
cd $dirday
mv ../../${dirday}/ft_*.LHZ.SAC .
if( `ls ../../${dirday} | wc -l` == 0 )rm -rf ../../${dirday}
ls ft_* > temp.ft
awk '{print "60 50 10 7.5 1 1  ", $1}' temp.ft > param.dat
/home/tianye/code/Programs/do_before_cor/filter4/filter4 param.dat 1
awk '{print "100 50 5 4 1 1  ", $1, $1"_10.0_50.0"}' temp.ft > param.dat
/home/tianye/code/Programs/do_before_cor/whiten_outphamp_C/whiten_rej_phamp_no_fil  param.dat /home/tianye/code/Programs/do_before_cor/whiten_outphamp_C/Stack.smooth
rm -f ft_*.SAC ft_*.SAC_10.0_50.0 param.dat temp.ft
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

mkdir -p COR

/home/tianye/code/Programs/CORR_REC/justCOR 3000 1

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
