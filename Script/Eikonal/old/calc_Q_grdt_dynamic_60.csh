#!/bin/csh
if ( $#argv != 1)then
echo "USAGE: get_event [SAC_path]"
exit
endif

cd /home/tianye/code/Programs/Attenuation/Alpha_from_ampg_btag_tvtl_test_coef/
cp driver_dynamic.C driver.C
make clean
make
set mod_pa='/home/tianye/Model'

cd $argv[1]

#foreach per ( 10 15 20 25 30 40 50 60 70 80 )
foreach per (60)
@ dis = 4 * $per
cd $per'sec_6snr_'$dis'dis'
rm -f ev_files.lst
foreach evinfo (`awk '{print $1"@"$2"@"$3}' ../event.loc`)
set ev=`echo $evinfo | cut -d@ -f1`
if( ! -e $ev'.ph.txt' ) continue
set lon=`echo $evinfo | cut -d@ -f2`
set lat=`echo $evinfo | cut -d@ -f3`
#set dist=`~/code/Programs/DIST/get_dist 41 250 $lat $lon d`
#if( `echo $dist | awk '{if($1<2000)print 1; else print 0}'` ) continue
#if( ! -e $ev'_am.txt' )then
#awk '{print $1,$2,$5}' $ev'.ph.txt' > $ev'_am.txt'
#endif
#/home/tianye/code/Programs/Earthquake/correct_2pi_v1 $ev'.ph.txt' $per 50
#if( ! -e $ev'.ph.txt_v1' )continue
#awk '{print $1,$2,$5}' $ev'.ph.txt_v1' > $ev'_am.txt_v1'
#cp $ev'.ph.txt' $ev'.ph.txt_v1'
#cp $ev'_am.txt' $ev'_am.txt_v1'
#C_plot_travel $ev'.ph.txt_v1' ~/region_TA 0.1
#C_plot_travel $ev'_am.txt_v1' ~/region_TA 0.1
#/home/tianye/code/Programs/Earthquake/correct_travel_time_curvature_v1_270 $ev $per 50 0.02 0.05
#rm -f $ev'.ph.txt_v1.HD' $ev'_am.txt_v1.HD'
#awk '{print $1,$2,$5}' $ev'.ph.txt_v2' > $ev'_am.txt_v2'
echo $ev'.ph.txt_v2' $ev'_am.txt_v2' $lon $lat >> ev_files.lst
end
/home/tianye/code/Programs/Attenuation/Alpha_from_ampg_btag_tvtl_test_coef/alpha_from_ampg_btag_tvtl ev_files.lst ../station.lst $mod_pa'/beta/beta_'$per
mkdir -p alpha_sm_0_result_-0.003_0.003_dynamic_2.0
mv 'alpha_'*'_map' alpha_sm_0_result_-0.003_0.003_dynamic_2.0
mv site_coef_map alpha_sm_0_result_-0.003_0.003_dynamic_2.0
cd alpha_sm_0_result_-0.003_0.003_dynamic_2.0
foreach file (alpha_raw_map alpha_prp_crd_map alpha_sit_crd_map)
awk '$4>30' $file > $file'_1'
~/code/Programs/mini_tools/Gauss_Smoothing_map $file'_1' 100
C_plot_travel $file'_1_smd' ~/region_TA 0.1
C_plot_input_region_res_csta '../alpha_'$per'.cpt' $file'_1_smd.HD' ~/region_TA 0.1 0_0
end
set file='alpha_crd_map'
awk '$4>30' $file > $file'_1'
foreach sm (100 300 500 800)
~/code/Programs/mini_tools/Gauss_Smoothing_map $file'_1' $sm
mv $file'_1_smd' $file'_1_smd_'$sm
C_plot_travel $file'_1_smd_'$sm ~/region_TA 0.1
C_plot_input_region_res_csta '../alpha_'$per'.cpt' $file'_1_smd_'$sm'.HD' ~/region_TA 0.1 0_0
end
foreach file (alpha_site_map alpha_prp_map)
awk '$4>30' $file > $file'_1'
C_plot_travel $file'_1' ~/region_TA 0.1
C_plot_input_region_res_csta '../alpha_'$per'.cpt' $file'_1.HD' ~/region_TA 0.1 0_0
end
cd ../..
end

