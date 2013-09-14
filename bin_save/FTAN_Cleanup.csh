#!/bin/csh
if ( $#argv != 1)then
echo "USAGE: "$0" [station dir]"
exit
endif

set sta = $argv[1]
rm -f $sta'/param_R.dat' $sta'/sac.lst'
### delete all .SAC_pos .SAC_neg .SAC_sym _?_AMP _?_DISP.? _amp_snr files
find $sta \( -name 'COR_'$sta'_*.SAC_???' -or -name 'COR_'$sta'_*.SAC*_?_AMP' -or -name 'COR_'$sta'_*.SAC*_?_DISP.?' -or -name 'COR_'$sta'_*.SAC*_amp_snr' \) -type f -print | xargs rm -f

