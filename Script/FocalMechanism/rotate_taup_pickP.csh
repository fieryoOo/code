#!/bin/csh
if ($#argv != 4) then
  echo "USAGE: rotate_taup_pickP.csh [sac_path] [ev.loc] [sta_lst] [arrival_file]"
  exit 1
endif

set ev=`awk '{print $1}' $argv[2]`
set elon=`awk '{print $2}' $argv[2]`
set elat=`awk '{print $3}' $argv[2]`
set ar_file=$argv[4]
set dep=8.
set ep_t0=1203603362.72700
cd $argv[1]

set sacf = ( " " " " " " " " " " )
rm -f Wells_file_dis_azi_itak_iinc_t0_tvt_pha_pol.txt
foreach info (`awk '{print $1"@"$2"@"$3}' $argv[3]`)

   set sta=`echo $info | cut -d@ -f1`
   set slon=`echo $info | cut -d@ -f2`
   set slat=`echo $info | cut -d@ -f3`

   set sacf[1]=`ls $ev*$sta*HE.sac | awk 'NR==1'`
   if ( `ls "$sacf[1]" | wc -l` == 0 ) continue
   set sacf[2]=`ls $ev*$sta*HN.sac | awk 'NR==1'`
   if ( `ls "$sacf[2]" | wc -l` == 0 ) continue
   set sacf[3]=`ls $ev*$sta*HZ.sac | awk 'NR==1'`
   if ( `ls "$sacf[3]" | wc -l` == 0 ) continue

   set sacf[4]=`echo $sacf[1] | sed s/'HE.sac'/'HR.sac'/g`
   set sacf[5]=`echo $sacf[1] | sed s/'HE.sac'/'HT.sac'/g`

/home/tianye/Software/sac/bin/sac << END
r $sacf[1] $sacf[2]
rotate to GCP
w $sacf[4] $sacf[5]
quit
END

   set dis=`/home/tianye/code/Programs/DIST/get_dist $elat $elon $slat $slon d`
   set azi=`/home/tianye/code/Programs/DIST/get_dist $elat $elon $slat $slon a1`

   /home/tianye/Software/TauP/TauP-2.1.1/bin/taup_time -mod prem -ph P,Pn,Pg -h $dep -km $dis | awk 'NR>5 && NF>8' > temp1
   set flag=0
   foreach tmp_info (`awk '{print $3"_"$4"_"$6"_"$7}' temp1`)
      set ph=`echo $tmp_info | cut -d_ -f1`
      if ( $ph != "P" ) continue
      set t0=`echo $tmp_info | cut -d_ -f2`
      set itak=`echo $tmp_info | cut -d_ -f3`
      set iinc=`echo $tmp_info | cut -d_ -f4`
      if ( `awk -v t0=$t0 -v itak=$itak '$3!="P" && ($4-t0)**2<0.01**2 && ($6-itak)**2<0.01**2' temp1 | wc -l` ) continue
      set flag=1
      break
   end
   if( $flag == 0 ) continue

   set ar_info=`awk -v sta=$sta '$1==sta{print $2"_"$7"_"$8"_"$21}' $ar_file`
   set tvt=`echo $ar_info | cut -d_ -f1 | awk -v t0=$ep_t0 '{print $1-t0}'`
   set ch=`echo $ar_info | cut -d_ -f2`
   set pha=`echo $ar_info | cut -d_ -f3`
   set pol=`echo $ar_info | cut -d_ -f4`
   set f_name=`echo $sacf[1] | sed s/'BHE.sac'/$ch'.sac'/g`

   echo $f_name $dis $azi $itak $iinc $t0 $tvt $pha $pol >> Wells_file_dis_azi_itak_iinc_t0_tvt_pha_pol.txt
   
end
