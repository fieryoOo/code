rm -f Wells_file_dis_azi_itak_tvt_pol_R_amp.txt
rm -f Plot_Focal_amp.txt
set MT = (1.46e+24 -3.14e+24 6.73e+24 6.10e+23 -1.31e+24 -8.19e+24)
foreach info (`awk '{print $3"@"$4"@"$8"@"$9"@"$1"@"$2"@"$7}' Wells_file_dis_azi_itak_iinc_t0_tvt_pha_pol.txt`)
set PH=`echo $info | cut -d@ -f3`
if( $PH != "P" ) continue
set pol=`echo $info | cut -d@ -f4 | cut -c1`
if( $pol == "c" ) then
set R=1
else if ( $pol == "d" ) then
set R=-1
else continue
endif
set azi=`echo $info | cut -d@ -f1`
set itak=`echo $info | cut -d@ -f2`
set file=`echo $info | cut -d@ -f5`
set dis=`echo $info | cut -d@ -f6`
set tvt=`echo $info | cut -d@ -f7`
set gxx=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print (sin($1*pi/180.))**2*(cos($2*pi/180.))**2}'`
set gxy=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print (sin($1*pi/180.))**2*sin(2*$2*pi/180.)}'`
set gyy=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print (sin($1*pi/180.))**2*(sin($2*pi/180.))**2}'`
set gxz=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print sin(2*$1*pi/180.)*cos($2*pi/180.)}'`
set gyz=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print sin(2*$1*pi/180.)*sin($2*pi/180.)}'`
set gzz=`echo $itak | awk 'BEGIN{pi=3.141592654}{print (cos($1*pi/180.))**2}'`

#echo $gxx $gxy $gyy $gxz $gyz $gzz | awk -v R=$R '{print -$1*R,-$2*R,-$3*R,-$4*R,-$5*R,-$6*R,$1*R,$2*R,$3*R,$4*R,$5*R,$6*R,-1e-5}' >> temp1
set amp=`echo $gxx $gxy $gyy $gxz $gyz $gzz $MT[1] $MT[2] $MT[3] $MT[4] $MT[5] $MT[6] | awk '{print $1*$7+$2*$8+$3*$9+$4*$10+$5*$11+$6*$12}'` 
echo $file $dis $azi $itak $tvt $pol $R $amp >> Wells_file_dis_azi_itak_tvt_pol_R_amp.txt
echo $azi $itak $amp | awk 'BEGIN{pi=3.141592654}{print sin($2*pi/180.)*sin($1*pi/180.), sin($2*pi/180.)*cos($1*pi/180.), $3}' >> Plot_Focal_amp.txt
end

