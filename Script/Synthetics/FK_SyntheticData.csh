### set up source geometry
set clon = 0.
set clat = 0.
set a = 30.
set b = 30.
set ang = -30.
#echo $clon $clat $a $b $ang | awk '{clon=$1; clat=$2; a=$3; b=$4; ang=$5; pi=3.1415926536; for(i=0;i<360;i++){theta = i*pi/180.; r = a*b/((b*cos(theta))**2+(a*sin(theta))**2)**0.5; print r*cos(theta+ang*pi/180.)+clon, r*sin(theta+ang*pi/180.)+clat, sin(theta/2.)*1e18, 0.05 } }' > source.lst
echo $clon $clat $a $b $ang | awk '{clon=$1; clat=$2; a=$3; b=$4; ang=$5; pi=3.1415926536; for(i=0;i<360;i++){theta = i*pi/180.; r = a*b/((b*cos(theta))**2+(a*sin(theta))**2)**0.5; print r*cos(theta+ang*pi/180.)+clon, r*sin(theta+ang*pi/180.)+clat, sin(theta/2.)*1e18, 0.02 } }' > source.lst
exit

set npts = 84001
set outf = seed_syn.lst
set year = 2000
set nmonth = 1
foreach month (JAN FEB MAR APR MAY JUN)
   set dir = $year'.'$month
   mkdir -p $dir
   set nd = 31
   if( $month == "APR" || $month == "JUN" || $month == "SEP" || $month == "NOV" ) then
      set nd = 30
   else if ( $month == "FEB" ) then
      set nd = 29
   endif
   set day = 1
   while( $day <= $nd )
      mkdir -p $dir"/"$dir"."$day
      echo $dir"."$day
      echo $dir $day | awk '{for(ista=0;ista<6;ista++)print ista*2, 0, $1"/"$1"."$2"/"$1"."$2".STA"ista}' > sta_lst
      echo $dir"."$day".seed" $year $nmonth $day >> $outf
      foreach file (`echo $dir $day | awk '{for(ista=0;ista<6;ista++)print $1"/"$1"."$2"/ft_"$1"."$2".STA"ista".BHZ.SAC_rec"}'`)
         echo 0 $npts > $file
      end
      @ day ++
      NoiseSynthetic source.lst sta_lst
   end #while day
   @ nmonth ++
end #foreach month

csh filt.csh
