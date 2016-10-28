#!/bin/bash
shopt -s expand_aliases
source ~/.bashrc

Initialize() {
	fhlen=0.008
   fcleanup=filetobedel.txt
   ftable=FTAN_files_fund.lst
   #ModDir=/mtera/tianye/ASN_OBS_BHZ/PaperResults/DISP_Pred/SV_allpath_model
	ModDir=/work2/tianye/ASN_OBS_BHZ/PaperResults/DISP_Pred/SV_allpath_model
	lFTAN=1_DISP.0	# 2_DISP.1 (after phase match filter)
	iFTAN=`echo $lFTAN | cut -d_ -f1`
	if [ ! -e ${ModDir} ]; then
		echo "Invalid mod path: "${ModDir}
		exit
	fi
   FTANexe='/home/tianye/code/Programs/FTAN_amp_snr_predgrv/aftani_c_pgl_amp'
	if [ ! -e ${FTANexe} ]; then
		echo "Invalid FTANexe path: "${FTANexe}
		exit
	fi
	#&piover4,&vmin,&vmax,&tmin,&tmax,&tresh,&ffact,&taperl,&snr,&fmatch, &tresh2,&ffact2,&taperl2,&snr2,&fmatch2,name,&flag
   param='-1 0.15 4.5 0.8 35. 30. 2.0 3. 0.99 8. 15. 1.5 3.0 0.6 3.'
   \rm -f $ftable
}

GetNames() {
   ls ${sacf}_* 2>/dev/null | grep -v '\.ps' | xargs rm -f
   sacn=$sacf'_neg'
   sacp=$sacf'_pos'
   echo $sacp $sacn >> $fcleanup
   sac=`echo $sacf | awk -F/ '{print $NF}'`
   sta2=`echo $sac | cut -d_ -f3 | cut -d. -f1`
}

CheckPathType_old() {
   locs=`awk -v sta1=$sta1 -v sta2=$sta2 '{if($1==sta1){lon1=$2; lat1=$3; if(lon1<0){lon1+=360.}}else if($1==sta2){lon2=$2; lat2=$3; if(lon2<0.){lon2+=360.} print lon1,lat1,lon2,lat2; exit}}' $fsta`
   ### decide which model to be used according to the path location
   ### for any path with at least one station in 233.5<lon<235.7, use neither model
   ### assumes good oceanic model near the ridge for paths with both station lon<233.5 and paths with over 70% length in the ocean
   ### assumes averaged continental model would be good enough for any continental paths and pths with over 71% length on the continent
   pathtype=`echo $locs | awk '{if($1<$3){lon1=$1;lon2=$3;}else{lon1=$3;lon2=$1} bdl=233.5; bdr=235.7; fctr=(bdr-lon1)/(lon2-lon1); if(lon1<bdr&&lon1>bdl || lon2<bdr&&lon2>bdl){print 2}else if(fctr>0.7){print 0}else if(fctr<0.3){print 1}else{print 2}}'`
   if [ $pathtype == 0 ]; then echo "   oceanic path";
   elif [ $pathtype == 1 ]; then echo "   continental path";
   else echo "   combined path"
   fi
}

CheckPathType() {
	s2type=`grep $sta2 $fsta | awk '{print $4}'`
	pathtype=${s1type}${s2type}
	if [ $s1type == "D" ] || ( [ $s1type == "S" ] && [ $s2type == "C" ] ); then
		pathtype=${s2type}${s1type}
	fi
	echo "path_type = "$pathtype
}

GetPosNegSAC() {
   if [ `saclst DEPMIN f $sacf | awk '{if($2=="nan"){print 0}else{print 1}}'` == 0 ]; then
      echo "   empty SAC file!"
      continue
   fi
sac << END > /dev/null
cut -3000 0
r $sacf
reverse
ch b 0
w $sacn
cut 0 3000
r $sacf
w $sacp
quit
END
echo "   GetPosNegSAC is done.."
}

RunFTAN_old() {
   if [ $pathtype == 0 ]; then
      DispGrv=${ModDir}/${sta1}_${sta2}.R.grv
      DispPhv=${ModDir}/${sta1}_${sta2}.R.phv
   elif [ $pathtype == 1 ]; then
      DispGrv=${ModDir}/Continental.R.grv
      DispPhv=${ModDir}/Continental.R.phv
   elif [ $pathtype == 2 ]; then
      DispGrv='I_Dont_believe_this_name_could_exist'
      DispPhv=${ModDir}/J42A_J50A.R.phv
   else
      echo "Unknow path type: "$pathtype
      exit
   fi
   if [ $pathtype == 0 ] && ( [ ! -e $DispGrv ] || [ ! -e $DispPhv ] ); then
      echo "   No predictions found for path "${sta1}_${sta2}". Using J42A_J50A instead"
      DispGrv=${ModDir}/J42A_J50A.R.grv
      DispPhv=${ModDir}/J42A_J50A.R.phv
   fi
   echo -e $sacp"\n"$sacn | awk -v param="$param" '{print param,$1,"0" }' > param_R.dat
   echo $FTANexe param_R.dat $DispPhv $DispGrv ${fhlen} ${iFTAN}
   $FTANexe param_R.dat $DispPhv $DispGrv ${fhlen} ${iFTAN} >& /dev/null
   if [ $? != 0 ]; then
      echo "   ERROR: FTAN failed: "$FTANexe" param_R.dat "$DispPhv" "$DispGrv" "${fhlen}" "${iFTAN}
      exit
   else
      echo "   FTAN is done.."
   fi
}

RunFTAN() {
	# CC CS CD SS SD DD
	if [ $pathtype == "CC" ]; then
		prename=Cont_Cont.R
	elif [ $pathtype == "CS" ]; then
		prename=Cont_Shal.R
	elif [ $pathtype == "CD" ]; then
		prename=Cont_Deep.R
	else # oceannic, check fdisp_pred for the path first
		prename=${sta1}_${sta2}.R
		if [ ! -e ${ModDir}/${prename}.grv ]; then
			prename=${sta2}_${sta1}.R
		fi
		if [ ! -e ${ModDir}/${prename}.grv ]; then # not found, use the preds for the pathtype instead
			echo -n "Warning(RunFTAN): path prediction "${ModDir}"/"${prename}".grv not found. "
			if [ $pathtype == "SS" ] || [ $pathtype == "SD" ]; then
				prename=Shal_Deep.R
			else
				prename=Deep_Deep.R
			fi
			echo "Use "${ModDir}"/"${prename}".grv instead!"
		fi
	fi

	DispGrv=${ModDir}/${prename}.grv
	DispPhv=${ModDir}/${prename}.phv

   if [ ! -e $DispGrv ] || [ ! -e $DispPhv ]; then
      echo "Error(RunFTAN): predictions not exist: "${DispGrv}" "${DispPhv}
		exit
   fi

   echo -e $sacp"\n"$sacn | awk -v param="$param" '{print param,$1,"0" }' > param_R.dat
   echo $FTANexe param_R.dat $DispPhv $DispGrv ${fhlen} ${iFTAN}
   $FTANexe param_R.dat $DispPhv $DispGrv ${fhlen} ${iFTAN} >& /dev/null
   if [ $? != 0 ]; then
      echo "   ERROR: FTAN failed: "$FTANexe" param_R.dat "$DispPhv" "$DispGrv" "${fhlen}" "${iFTAN}
      exit
   else
      echo "   FTAN is done.."
   fi
}

ListResults() {
	#echo "debug in ListResults: "${ftable}" "${sacp}_${lFTAN}" "${sacp}_amp_snr" "${sacn}_${lFTAN}" "${sacn}_amp_snr
   if [ -e ${sacp}_${lFTAN} ] && [ -e ${sacp}_amp_snr ] && [ -e ${sacn}_${lFTAN} ] && [ -e ${sacn}_amp_snr ]; then
      dnum=`saclst user0 f $sacf | awk '{print $2}'`
      echo $sta1 $sta2 ${sacp}_${lFTAN} ${sacp}_amp_snr ${sacn}_${lFTAN} ${sacn}_amp_snr $dnum $DispGrv $DispPhv >> $ftable
   fi
}

CleanUp() {
   if [ -e $fcleanup ]; then
      more $fcleanup | xargs -n 1 sh -c 'rm -f $0'
      \rm -f $fcleanup
   fi
}

PickPath() {
   if [ $sacf != ${1}/COR_${1}_${2}.SAC ]; then
      continue
   fi
}

### main starts here ###
fsta=station-type.lst
if [ ! -e $fsta ]; then
	echo $fsta" not found"
	exit
fi

Initialize
while read sta1 lon lat s1type stmp; do
	#if [ $sta1 != "G03A" ]; then continue; fi
   for sacf in `ls ${sta1}/COR_${sta1}_*.SAC 2>/dev/null`; do
      #PickPath G03D G05D
      GetNames
		#if [ $sta2 != "J23A" ]; then continue; fi
		if [ $sta1 == $sta2 ]; then continue; fi
      echo "Working on path "${sta1}_${sta2}...
      CheckPathType
      GetPosNegSAC
      RunFTAN
      ListResults # check FTAN results and write into file table
   done
done < $fsta
CleanUp

