#!/bin/csh
if ( $#argv != 1 ) then
   echo "Usage: "$0" [station.lst]"
   exit
endif

### station list ###
set stalst = $argv[1]
if( ! -e $stalst ) then
   echo $stalst" not found!"
   exit
endif
###
set table = '/mtera/tianye/ASN_OBS_BHZ/SAC/STACK/FTAN_files_fund.lst'
rm -f $table
foreach sta (`awk '{print $1}' $stalst`)
#foreach sta (M03A)
echo "Working on station "$sta"..."
   #foreach sacf ( `ls ${sta}/COR_${sta}_J39A.SAC` )
   foreach sacf (`ls ${sta}/COR_${sta}_*.SAC`)
      /mtera/tianye/ASN_OBS_BHZ/PaperResults/DISP_Pred/Do_FTAN_fund $sacf $stalst $table
   end #sacf
end #sta

~/Programs/Disp/Investigate_FTAN_quality $stalst $table 1.2 5.5 0.03 0.015
~/Programs/Disp/Investigate_FTAN_quality $stalst $table 6. 20. 0.03 0.015

rm -f per.lst
foreach per (2. 4. 6. 8. 10. 14. 18. 22.)
echo -n $per" " >> per.lst
end
/home/tianye/code/Programs/Disp/Investigate_FTAN_wpred $stalst $table per.lst 0.35 0.15
rm -f per.lst
