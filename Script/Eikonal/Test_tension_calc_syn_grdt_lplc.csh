#!/bin/csh
if ( $#argv != 1)then
echo "USAGE: calc_syn_grdt_lplc.csh [scale(half per) in deg]"
exit
endif

###parameters
set clon=-114.8758
set clat=41.1003
set per=`echo $argv[1] | awk '{print $1*2}'`
set noise=0.02
set ts=0.20
set bs=20
###extract and shift station locations
mkdir -p per$per
cd per$per
awk -v clon=$clon -v clat=$clat '($2-clon)**2<49. && ($3-clat)**2<49. {print $2-clon, $3-clat, $1}' ../station.lst > syn_sta.lst
###produce cpt files
cp ../f.cpt .
awk -v per=$per 'BEGIN{c=(2*3.141592654/per)**2}{a=$1;b=$5;if(NR<=6){a=a*c;b=b*c}print a,$2,$3,$4,b,$6,$7,$8}' ../grdt.cpt > grdt.cpt
awk -v per=$per 'BEGIN{c=(2*3.141592654/per)**2}{a=$1;b=$5;if(NR<=6){a=a*c;b=b*c}print a,$2,$3,$4,b,$6,$7,$8}' ../lplc.cpt > lplc.cpt
###theorectical values at stations
awk -v per=$per 'BEGIN{pi=3.141592654}{c=2*pi/per; A=c*sqrt($1**2+$2**2); print $1,$2,sin(A)/A}' syn_sta.lst > 'syn_f_'$per'_noise0.txt'
awk -v per=$per 'BEGIN{pi=3.141592654}{c=2*pi/per; A=c*sqrt($1**2+$2**2); print $1,$2,c*c*(cos(A)/A-sin(A)/A/A)}' syn_sta.lst > 'syn_grdt_'$per'.txt'
awk -v per=$per 'BEGIN{pi=3.141592654}{c=2*pi/per; A=c*sqrt($1**2+$2**2); print $1,$2,c*c*(-(A**2-1)*sin(A)-A*cos(A))/(A**3)}' syn_sta.lst > 'syn_lplc_'$per'.txt'
###interporlate for syn original wave with and without blockmean
/home/tianye/code/Programs/mini_tools/add_noise 'syn_f_'$per'_noise0.txt' 'syn_f_'$per'_noise'$noise'.txt' $noise
foreach file ( 'syn_f_'$per'_noise0.txt' 'syn_f_'$per'_noise'$noise'.txt' )
../C_plot_travel $file ../region 0.1 0.00 $bs
../C_plot_travel $file ../region 0.1 0.00 1
end
###interporlate for syn gradient with and without blockmean
foreach file ( 'syn_grdt_'$per'.txt' )
../C_plot_travel $file ../region 0.1 $ts $bs
end
###interporlate for syn laplacian no tension with blockmean
foreach file ( 'syn_lplc_'$per'.txt' )
../C_plot_travel $file ../region 0.1 0.00 $bs
end
###theorectical HD surface
set HD_f = 'syn_grdt_'$per'.txt_ts'$ts'_bs'$bs'.HD'
awk -v per=$per 'BEGIN{pi=3.141592654}{c=2*pi/per; A=c*sqrt($1**2+$2**2); if(A==0){B=1}else{B=sin(A)/A} print $1,$2,B}' $HD_f > 'theo_f_'$per'.txt.HD'
awk -v per=$per 'BEGIN{pi=3.141592654}{c=2*pi/per; A=c*sqrt($1**2+$2**2); if(A==0){B=0}else{B=c*c*(cos(A)/A-sin(A)/A/A)} print $1,$2,B}' $HD_f > 'theo_grdt_'$per'.txt.HD'
awk -v per=$per 'BEGIN{pi=3.141592654}{c=2*pi/per; A=c*sqrt($1**2+$2**2); if(A==0){B=-2.*c*c/3.}else{B=c*c*(-(A**2-1)*sin(A)-A*cos(A))/(A**3)} print $1,$2,B}' $HD_f > 'theo_lplc_'$per'.txt.HD'
###plot theo and syn maps
echo 'theo_f_'$per'.txt.HD' 'f.cpt' > temp.lst
echo 'syn_f_'$per'_noise0.txt_ts0.00_bs'$bs'.HD' 'f.cpt' >> temp.lst
echo 'theo_grdt_'$per'.txt.HD' 'grdt.cpt' >> temp.lst
echo 'syn_grdt_'$per'.txt_ts'$ts'_bs'$bs'.HD' 'grdt.cpt' >> temp.lst
echo 'theo_lplc_'$per'.txt.HD' 'lplc.cpt' >> temp.lst
echo 'syn_lplc_'$per'.txt_ts0.00_bs'$bs'.HD' 'lplc.cpt' >> temp.lst
../Plot_6.csh temp.lst 'theo_syn_'$per'.ps' ../region
###calculate gradient and laplacian with and without ts&bs for wavefields with and without noise 
foreach file ( 'syn_f_'$per'_noise0.txt' 'syn_f_'$per'_noise'$noise'.txt' )
../calc_grdt_lplc $file 0.00 1
../calc_grdt_lplc $file $ts $bs
end
###interporlate for calculated gradient with and without tension/noise
foreach file ( 'syn_f_'$per'_noise'$noise'.txt' 'syn_f_'$per'_noise0.txt' )
../C_plot_travel_nonm $file'_ts0.00_bs1_grdt' ../region 0.1 0.00 1
../C_plot_travel_nonm $file'_ts0.00_bs'$bs'_grdt' ../region 0.1 $ts $bs
end
###interporlate for calculated laplacian no tension with blockmean
foreach file ( 'syn_f_'$per'_noise'$noise'.txt' 'syn_f_'$per'_noise0.txt' )
#../C_plot_travel_nonm $file'_ts0.00_bs1_lplc' ../region 0.1 0.00 $bs
#../C_plot_travel_nonm $file'_ts'$ts'_bs'$bs'_lplc' ../region 0.1 0.00 $bs
C_plot_travel_nn $file'_ts0.00_bs1_lplc' ../region 0.1
C_plot_travel_nn $file'_ts'$ts'_bs'$bs'_lplc' ../region 0.1
end
###plot maps computed without noise
echo 'syn_f_'$per'_noise0.txt_ts0.00_bs1.HD' 'f.cpt' > temp.lst
echo 'syn_f_'$per'_noise0.txt_ts0.00_bs'$bs'.HD' 'f.cpt' >> temp.lst
echo 'syn_f_'$per'_noise0.txt_ts0.00_bs1_grdt.HD' 'grdt.cpt' >> temp.lst
echo 'syn_f_'$per'_noise0.txt_ts0.00_bs'$bs'_grdt.HD' 'grdt.cpt' >> temp.lst
echo 'syn_f_'$per'_noise0.txt_ts0.00_bs1_lplc.HD' 'lplc.cpt' >> temp.lst
echo 'syn_f_'$per'_noise0.txt_ts'$ts'_bs'$bs'_lplc.HD' 'lplc.cpt' >> temp.lst
../Plot_6_2.csh temp.lst 'noise0_tsbs_'$per'.ps' ../region
###plot maps computed with noise
echo 'syn_f_'$per'_noise'$noise'.txt_ts0.00_bs1.HD' 'f.cpt' > temp.lst
echo 'syn_f_'$per'_noise'$noise'.txt_ts0.00_bs'$bs'.HD' 'f.cpt' >> temp.lst
echo 'syn_f_'$per'_noise'$noise'.txt_ts0.00_bs1_grdt.HD' 'grdt.cpt' >> temp.lst
echo 'syn_f_'$per'_noise'$noise'.txt_ts0.00_bs'$bs'_grdt.HD' 'grdt.cpt' >> temp.lst
echo 'syn_f_'$per'_noise'$noise'.txt_ts0.00_bs1_lplc.HD' 'lplc.cpt' >> temp.lst
echo 'syn_f_'$per'_noise'$noise'.txt_ts'$ts'_bs'$bs'_lplc.HD' 'lplc.cpt' >> temp.lst
../Plot_6_2.csh temp.lst 'noise'$noise'_tsbs_'$per'.ps' ../region

