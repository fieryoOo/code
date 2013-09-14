#!/bin/csh
if ( $#argv != 3)then
echo "USAGE: calc_agrdt_tlplc.csh [in_ph_file] [per] [in_beta_file]"
exit
endif

set file=$argv[1]
set ev=`echo $file | cut -d. -f1`
set per=$argv[2]
set clon=`head -1 $file | awk '{print $1}'`
set clat=`head -1 $file | awk '{print $2}'`

/home/tianye/code/Programs/Earthquake/correct_2pi_v1 $file $per 50
if( ! -e $file'_v1' )then
echo "no enough station points after 2pi correction"
exit
endif

awk '{print $1,$2,$5}' $file'_v1' > $ev'_am.txt_v1'
C_plot_travel $ev'.ph.txt_v1' ~/region_TA 0.1
C_plot_travel $ev'_am.txt_v1' ~/region_TA 0.1
/home/tianye/code/Programs/Earthquake/correct_travel_time_curvature_v1_270 $ev $per 50 0.02 0.05
awk '{print $1,$2,$5}' $file'_v2' > $ev'_am.txt_v2'

/home/tianye/code/Programs/Attenuation/csta_atn_aplf_term_earthquake $file'_v2' $ev'_am.txt_v2' $argv[3] $clon $clat

TXT2CPT_mid0 $ev'_alpha'
foreach file ( $ev'_alpha' $ev'_alpha_decay' $ev'_alpha_prp' $ev'_alpha_site' )
C_plot_travel $file ~/region_TA 0.1
C_plot_input_region_res_csta '../'$ev'_alpha.cpt' $file'.HD' ~/region_TA 0.1 0_0
end

