#!/bin/csh
if ( $#argv != 3)then
echo "USAGE: get_event [SAC_path] [event.info] [station.lst]"
exit
endif

set sm_hdis = 0
set p_coef = 1.0
set s_coef = 1.4

cd /home/tianye/code/Programs/Attenuation/Alpha_from_ampg_btag_tvtl_cap/
cp driver_40_dynamic.C driver.C
make clean
make
set mod_pa='/home/tianye/Model'

cd $argv[1]

foreach per (40)
@ dis = 4 * $per
cd $per'sec_6snr_'$dis'dis'
rm -f ev_files.lst
foreach evinfo (`awk '{print $1"@"$2"@"$3"@"$5"@"$6}' $argv[2]`)
set ev=`echo $evinfo | cut -d@ -f1`
if( ! -e $ev'.ph.txt' ) continue
set lon=`echo $evinfo | cut -d@ -f2`
set lat=`echo $evinfo | cut -d@ -f3`
set dep=`echo $evinfo | cut -d@ -f4`
set mag=`echo $evinfo | cut -d@ -f5`
#set dist=`~/code/Programs/DIST/get_dist 41 250 $lat $lon d`
#if( `echo $dist | awk '{if($1>5000)print 1; else print 0}'` ) continue
if( ! -e $ev'_am.txt' )then
awk '{print $1,$2,$5}' $ev'.ph.txt' > $ev'_am.txt'
endif
rm -f $ev'.ph.txt_v1' $ev'.ph.txt_v2' $ev'_am.txt_v1' $ev'_am.txt_v2'
/home/tianye/code/Programs/Earthquake/correct_2pi_v1 $ev'.ph.txt' $per 50
if( ! -e $ev'.ph.txt_v1' )continue
awk '{print $1,$2,$5}' $ev'.ph.txt_v1' > $ev'_am.txt_v1'
C_plot_travel $ev'.ph.txt_v1' ~/region_TA 0.1
C_plot_travel $ev'_am.txt_v1' ~/region_TA 0.1
/home/tianye/code/Programs/Earthquake/correct_travel_time_curvature_v1_270 $ev $per 50 0.02 0.05
rm -f $ev'.ph.txt_v1.HD' $ev'_am.txt_v1.HD'
if( ! -e $ev'.ph.txt_v2' )continue
awk '{print $1,$2,$5}' $ev'.ph.txt_v2' > $ev'_am.txt_v2'
echo $ev'.ph.txt_v2' $ev'_am.txt_v2' $lon $lat $dep $mag >> ev_files.lst
end

/home/tianye/code/Programs/Attenuation/Alpha_from_ampg_btag_tvtl_cap/alpha_from_ampg_btag_tvtl ev_files.lst $argv[3] $mod_pa'/beta/beta_'$per $sm_hdis $p_coef $s_coef

mkdir -p 'alpha_sm_'$sm_hdis'_result_'$p_coef'_'$s_coef
mv 'alpha_'*'_map' 'alpha_sm_'$sm_hdis'_result_'$p_coef'_'$s_coef
mv site_coef_map 'alpha_sm_'$sm_hdis'_result_'$p_coef'_'$s_coef
mv vel_map 'alpha_sm_'$sm_hdis'_result_'$p_coef'_'$s_coef
cd 'alpha_sm_'$sm_hdis'_result_'$p_coef'_'$s_coef
foreach file (alpha_crd_map alpha_raw_map alpha_prp_crd_map alpha_sit_crd_map)
#awk '$4>0 && $5>200' $file > $file'_1'
awk '$4>0 && $5>50' $file > $file'_1'
C_plot_travel_nn $file'_1' ~/region_TA 0.1
C_plot_input_region_res_csta '../alpha_'$per'.cpt' $file'_1.HD' ~/region_TA 0.1 0_0
~/code/Programs/mini_tools/Gauss_Smoothing_map_avg $file'_1' 100
mv $file'_1_smd' $file'_1_smd_100'
C_plot_travel_nn $file'_1_smd_100' ~/region_TA 0.1
C_plot_input_region_res_csta '../alpha_'$per'.cpt' $file'_1_smd_100.HD' ~/region_TA 0.1 0_0
end
set file='alpha_crd_map'
awk '$4>0 && $5>50' $file > $file'_1'
foreach sm ( 200 300 )
~/code/Programs/mini_tools/Gauss_Smoothing_map_avg $file'_1' $sm
mv $file'_1_smd' $file'_1_smd_'$sm
C_plot_travel_nn $file'_1_smd_'$sm ~/region_TA 0.1
C_plot_input_region_res_csta '../alpha_'$per'.cpt' $file'_1_smd_'$sm'.HD' ~/region_TA 0.1 0_0
end
foreach file (alpha_site_map alpha_prp_map)
awk '$4>0 && $5>50' $file > $file'_1'
C_plot_travel_nn $file'_1' ~/region_TA 0.1
C_plot_input_region_res_csta '../alpha_'$per'.cpt' $file'_1.HD' ~/region_TA 0.1 0_0
end
cd ../..
end

