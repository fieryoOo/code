#!/bin/csh
if ( $#argv != 4)then
echo "usage: Compare_trvt_amp [sac_path_1] [sac_path_2] [longitude] [latitude]"
exit
endif

set net1=`echo $argv[1] | awk -F/ '{print $NF}'`
if ( ! $#net1 ) set net1=`echo $argv[1] | awk -F/ '{print $(NF-1)}'`
set net2=`echo $argv[2] | awk -F/ '{print $NF}'`
if ( ! $#net2 ) set net2=`echo $argv[2] | awk -F/ '{print $(NF-1)}'`
#echo $net1,$net2
#exit

foreach period ( 30 50 80 )
cd $argv[1]
echo 'working on '$net1' for '$period'sec...'
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis'
rm -f 'trvt_amp_term_'$period's_'$argv[3]'_'$argv[4]
foreach event_am (`ls 20*_am_laplace.txt.HD`)
set event_sl=`echo $event_am | cut -d_ -f1`
set event_sl=`echo 'slow_azi_'$event_sl'.txt.HD'`
set trvt_term=`awk -v long=$argv[3] -v lati=$argv[4] '(long-$1)^2+(lati-$2)^2<0.000001 {printf "%.5f",$3^2}' $event_sl`
set amp_term=`awk -v long=$argv[3] -v lati=$argv[4] '(long-$1)^2+(lati-$2)^2<0.000001 {printf "%.5f",$3}' $event_am`
set amp_temp=`echo "sqrt(($amp_term*100000)^2)*2/2" | bc`
if ( `echo "$trvt_term*100000" | bc` == 0 | $amp_temp == 0 )continue
set v1000=`echo $trvt_term | awk '{printf "%.0f",1000/sqrt($1)}'`
if ( $v1000 < 3000 | $v1000 > 4500 | $amp_temp > 1000 )continue
echo $trvt_term $amp_term >> 'trvt_amp_term_'$period's_'$argv[3]'_'$argv[4]
#exit
end

cd $argv[2]
echo 'working on '$net2' for '$period'sec...'
@ dis = $period * 12
cd $period'sec_10snr_'$dis'dis'
rm -f 'trvt_amp_term_'$period's_'$argv[3]'_'$argv[4]
foreach event_am (`ls 20*_am_laplace.txt.HD`)
set event_sl=`echo $event_am | cut -d_ -f1`
set event_sl=`echo 'slow_azi_'$event_sl'.txt.HD'`
set trvt_term=`awk -v long=$argv[3] -v lati=$argv[4] '(long-$1)^2+(lati-$2)^2<0.000001 {printf "%.5f",$3^2}' $event_sl`
set amp_term=`awk -v long=$argv[3] -v lati=$argv[4] '(long-$1)^2+(lati-$2)^2<0.000001 {printf "%.5f",$3}' $event_am`
set amp_temp=`echo "sqrt(($amp_term*100000)^2)*2/2" | bc`
if ( `echo "$trvt_term*100000" | bc` == 0 | $amp_temp == 0 )continue
set v1000=`echo $trvt_term | awk '{printf "%.0f",1000/sqrt($1)}'`
if ( $v1000 < 3000 | $v1000 > 4500 | $amp_temp > 1000 )continue
echo $trvt_term $amp_term >> 'trvt_amp_term_'$period's_'$argv[3]'_'$argv[4]
#exit
end
set file1=$argv[1]'/'$period'sec_10snr_'$dis'dis/trvt_amp_term_'$period's_'$argv[3]'_'$argv[4]
set file2=$argv[2]'/'$period'sec_10snr_'$dis'dis/trvt_amp_term_'$period's_'$argv[3]'_'$argv[4]
set avg=`awk 'BEGIN{a=0}{a+=$1}END{print a}' $file1`
set N1=`more $file1 | wc -l`
set avg=`awk -v a=$avg -v N1=$N1 '{a+=$1}END{printf "%.4f",a/(N1+NR)}' $file2`
gnuplot << END
f$net1(x)=a1*x+b1
fit f$net1(x) '$file1' using 2:1 via a1,b1
f$net2(x)=a2*x+b2
fit f$net2(x) '$file2' using 2:1 via a2,b2
set print 'temp.gnp'
print a1,' * x + ',b1
print a2,' * x + ',b2
set print
END
gnuplot << END
plot '$file1' using 2:1 title '$net1 $period sec $argv[3] $argv[4]',\
'$file2' using 2:1 title '$net2 $period sec $argv[3] $argv[4]',\
`awk '{if(NR==1){print}}' temp.gnp`,\
`awk '{if(NR==2){print}}' temp.gnp`
set xrange [-0.01:0.01]
set yrange [$avg-0.01:$avg+0.01]
set term postscript enhanced color
set output '$file1.ps'
replot
set output
set term x11
END
mv $file1.ps /home/tianye/data_Eikonal/results/trvt_amp_terms

end

