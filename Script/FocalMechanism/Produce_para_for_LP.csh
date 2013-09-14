set in_file='/media/WORK/tianye/Events_Wells/20080221141604/Wells_file_dis_azi_itak_tvt_pol_R_amp_57.txt'
rm -f temp1
foreach info (`awk 'NR<200{print $3"@"$4"@"$6}' $in_file`)
#set PH=`echo $info | cut -d@ -f3`
#if( $PH != "P" ) continue
set pol=`echo $info | cut -d@ -f3 | cut -c1`
if( $pol == "c" ) then
set R=1
else if ( $pol == "d" ) then
set R=-1
else continue
endif
set azi=`echo $info | cut -d@ -f1`
set itak=`echo $info | cut -d@ -f2`
set gxx=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print (sin($1*pi/180.))**2*(cos($2*pi/180.))**2}'`
set gxy=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print (sin($1*pi/180.))**2*sin(2*$2*pi/180.)}'`
set gyy=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print (sin($1*pi/180.))**2*(sin($2*pi/180.))**2}'`
set gxz=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print sin(2*$1*pi/180.)*cos($2*pi/180.)}'`
set gyz=`echo $itak $azi | awk 'BEGIN{pi=3.141592654}{print sin(2*$1*pi/180.)*sin($2*pi/180.)}'`
set gzz=`echo $itak | awk 'BEGIN{pi=3.141592654}{print (cos($1*pi/180.))**2}'`

echo $gxx $gxy $gyy $gxz $gyz $gzz | awk -v R=$R '{print -$1*R,-$2*R,-$3*R,-$4*R,-$5*R,-$6*R,$1*R,$2*R,$3*R,$4*R,$5*R,$6*R,-1e-10}' >> temp1
end
set N=`more temp1 | wc -l`
echo -1 12 $N > Wells_para_LP.txt
echo 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 >> Wells_para_LP.txt
cat temp1 >> Wells_para_LP.txt

