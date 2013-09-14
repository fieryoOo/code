#!/bin/csh
if ( $#argv != 4 )then
echo "usage: Plot_amp_dist_deg_map.csh [cent_sta] [sta_lst] [input_file] [deg_file as: min max]"
exit
endif

cd /home/tianye/data_center_sta/$argv[1]
set long=`grep -m1 $argv[1] $argv[2] | awk '{print $2}'`
set lati=`grep -m1 $argv[1] $argv[2] | awk '{print $3}'`

set cpt_txt=`echo $argv[3] | sed s/'neg'/'pos'/ | sed s/'sym'/'pos'/`
#set cpt_txt=`echo $argv[3] | sed s/'pos'/'neg'/ | sed s/'sym'/'neg'/`

mkdir -p amp_dist_$argv[4]
cp $argv[3] amp_dist_$argv[4]
cp $argv[3]'.HD' amp_dist_$argv[4]
cp $cpt_txt'.HD.cpt' amp_dist_$argv[4]
cp $argv[1]'.region' amp_dist_$argv[4]

cd amp_dist_$argv[4]
rm -f $argv[4]'.deg'
foreach deg (`more ../$argv[4]`)
echo $deg | awk '{print $1-180}' >> $argv[4]'.deg'
end

/home/tianye/code/Script/GMT/C_plot_kernel_input_angle $cpt_txt'.HD.cpt' $argv[3]'.HD' $argv[1]'.region' $long $lati $argv[4]'.deg'

foreach angle (`awk '{print $1"@"$2}' '../'$argv[4]`)
set deg_min=`echo $angle | cut -d@ -f1`
set deg_max=`echo $angle | cut -d@ -f2`
/home/tianye/code/Programs/CORRECT_ASN_AMP/with_days/amp_dist_angle $long $lati $argv[3] $deg_min $deg_max
#mv $argv[3]'_amp_dist_deg'$deg_min'_to_deg'$deg_max amp_dist_$argv[4]
end


rm -f temp.fit
foreach adfile (`ls $argv[3]'_amp_dist_deg'*"_to_deg"*`)
gnuplot << END
f(x)=a*x+b
a=-0.001
fit f(x) '$adfile' via a,b
set print 'temp.gnp'
print a,' * x + ',b
set print
END
more temp.gnp >> temp.fit
end

set a=`awk 'BEGIN{a=0}{a+=$1}END{print a/NR}' temp.fit`
set b=`awk 'BEGIN{a=0}{a+=$5}END{print a/NR}' temp.fit`
rm -f amp_dist_picked
rm -f amp_dist_discard
foreach adfile (`ls $argv[3]'_amp_dist_deg'*"_to_deg"*`)
awk -v a=$a -v b=$b '{if(sqrt((a*$1-$2+b)**2)/sqrt(a**2+1)<0.8){print}}' $adfile >> amp_dist_picked
awk -v a=$a -v b=$b '{if(sqrt((a*$1-$2+b)**2)/sqrt(a**2+1)>=0.8){print $1,$2,$3,sqrt((a*$1-$2+b)**2)/sqrt(a**2+1)}}' $adfile >> amp_dist_discard
end

rm -f temp.fit
#foreach adfile (`ls $argv[3]'_amp_dist_deg'*'_to_deg'*'_picked'`)
gnuplot << END
f(x)=a*x+b
a=-0.001
fit f(x) 'amp_dist_picked' via a,b
set print 'temp.fit'
print a,' * x + ',b
set print
END
#end

#rm -f temp.fit
#foreach adfile (`ls $argv[3]'_amp_dist_deg'*"_to_deg"*`)
#gnuplot << END
#f(x)=a*x+$b
#fit f(x) '$adfile' via a
#set print 'temp.gnp'
#print a,' * x + $b'
#set print
#END
#more temp.gnp >> temp.fit
#end
awk '{printf "%.5f %s %s %s %.2f\n",$1,$2,$3,$4,$5}' temp.fit > temp2.fit
mv temp2.fit temp.fit

#gnuplot << END
#plot '$adfile' title '`echo $argv[3] | sed s/'\_'/' '/g` $argv[4] $argv[5]',\
#`awk '{if(NR==1){print}}' temp.gnp`
#set term postscript enhanced color
#set output '$adfile.ps'
#replot
#set output
#set term x11
#END

#end
