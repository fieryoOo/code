#!/bin/csh
if ( $#argv != 2)then
echo "USAGE: plot_crct_sit_prp [site coef] [prp coef]"
exit
endif

set per=25
mkdir -p $argv[1]'sit_'$argv[2]'prp'
awk '{print $1"_"$2,$3}' alpha_raw_map_1 > $argv[1]'sit_'$argv[2]'prp/temp0'
awk '{print $1"_"$2,$3/1.5}' alpha_site_map_1 > $argv[1]'sit_'$argv[2]'prp/temps'
awk '{print $1"_"$2,$3}' alpha_prp_map_1 > $argv[1]'sit_'$argv[2]'prp/tempp'
cd $argv[1]'sit_'$argv[2]'prp'
join temp0 temps | awk -v coef=$argv[1] '{print $1,$2+coef*$3}' | awk -F_ '{print $1,$2}' > alpha_sit_crd_map
awk '{print $1"_"$2,$3}' alpha_sit_crd_map > tempsc
join temp0 tempp | awk -v coef=$argv[2] '{print $1,$2+coef*$3}' | awk -F_ '{print $1,$2}' > alpha_prp_crd_map
join tempsc tempp | awk -v coef=$argv[2] '{print $1,$2+coef*$3}' | awk -F_ '{print $1,$2}' > alpha_crd_map
#rm -f temp0 temps tempp tempsc

foreach file (alpha_crd_map alpha_prp_crd_map alpha_sit_crd_map)
~/code/Programs/mini_tools/Gauss_Smoothing $file 100
C_plot_travel $file'_smd' ~/region_TA 0.1
C_plot_input_region_res_csta '../../alpha_'$per'.cpt' $file'_smd.HD' ~/region_TA 0.1 0_0
end
