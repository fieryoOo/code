#!/bin/csh
if ( $#argv != 2)then
echo "USAGE: "$0" [SAC_path] [station.lst]"
exit
endif

set per = 14
cd $argv[1]
set dir = 'Timing_Error'
mkdir -p $dir
foreach station (`awk '{print $1}' $argv[2]`)
echo "Working on station "$station
rm -f $dir'/'$station'_misfit.txt'
foreach sacf (`ls $station/COR_J23A_*.SAC */COR_*_$station.SAC`)
set sta1 = `echo $sacf | cut -d/ -f1`
if( `echo $station $sta1 | awk '{if($1==$2)print 1}'` ) then
set sta2 = `echo $sacf | cut -d_ -f3 | cut -d. -f1`
set posf = $sacf"_pos_2_DISP.1"
set negf = $sacf"_neg_2_DISP.1"
else 
set sta2 = $sta1
set negf = $sacf"_pos_2_DISP.1"
set posf = $sacf"_neg_2_DISP.1"
endif
if( ! -e $posf ) continue
if( ! -e $negf ) continue
echo $station $sta2
set NRp=`awk '$3>'$per'{print $1}' $posf | head -1`
if( ! $NRp ) continue
if( $NRp == 0 ) continue
set NRn=`awk '$3>'$per'{print $1}' $negf | head -1`
if( ! $NRn ) continue
if( $NRn == 0 ) continue
set poss = $sacf"_pos_amp_snr"
set snrp = `awk 'NR>='$NRp' && NR<='$NRp+1'{printf "%f %f ", $1,$3}' $poss | awk '{print $2+('$per'-$1)/($3-$1)*($4-$2)}'`
if( `echo $snrp | awk '{if($1<6)print 1}'` ) continue
set negs = $sacf"_pos_amp_snr"
set snrn = `awk 'NR>='$NRn' && NR<='$NRn+1'{printf "%f %f ", $1,$3}' $negs | awk '{print $2+('$per'-$1)/($3-$1)*($4-$2)}'`
if( `echo $snrn | awk '{if($1<6)print 1}'` ) continue
set dist = `saclst dist f $sacf | awk '{print $2}'`
set velp = `awk 'NR>='$NRp' && NR<='$NRp+1'{printf "%f %f ", $3,$5}' $posf | awk '{print $2+('$per'-$1)/($3-$1)*($4-$2)}'`
set veln = `awk 'NR>='$NRn' && NR<='$NRn+1'{printf "%f %f ", $3,$5}' $negf | awk '{print $2+('$per'-$1)/($3-$1)*($4-$2)}'`
set misf = `echo $velp $veln | awk -v dis=$dist '{print dis/$1-dis/$2}'`
echo $misf $sta2 $velp $veln $dist >> $dir'/'$station'_misfit.txt'
end
end
