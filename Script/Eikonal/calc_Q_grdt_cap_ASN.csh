#!/bin/csh
if ( $#argv != 6)then
echo "USAGE: get_event [SAC_path] [c_sta.lst] [station.lst] [per.lst] [fix or dyn] [param_file]"
exit
endif

set mod_pa='/home/tianye/Model'

cd $argv[1]

foreach per (`more $argv[4]`)
set cpt_f=$argv[1]'/alpha_cpt/alpha_'$per'.cpt'
#@ dis = 4 * $per
#cd $per'sec_6snr_'$dis'dis'
cd 'Ph_Amp_Map_'$per'sec'
rm -f ev_files.lst
foreach evinfo (`awk '{print $1"@"$2"@"$3}' $argv[2]`)
set ev=`echo $evinfo | cut -d@ -f1`
if( ! -e $ev'_center_ph_amp_map' ) continue
set lon=`echo $evinfo | cut -d@ -f2`
set lat=`echo $evinfo | cut -d@ -f3`
#set dist=`~/code/Programs/DIST/get_dist 41 250 $lat $lon d`
#if( `echo $dist | awk '{if($1<1500)print 1; else print 0}'` ) continue
#if( ! -e $ev'_center_am_map' )then
#awk '{print $1,$2,$5,$6}' $ev'_center_ph_amp_map' > $ev'_center_am_map'
#endif # ! -e
#rm -f $ev'.ph.txt_v1' $ev'.ph.txt_v2' $ev'_am.txt_v1' $ev'_am.txt_v2'
#/home/tianye/code/Programs/Earthquake/correct_2pi_v1 $ev'_center_ph_amp_map' $per 50
#if( ! -e $ev'_center_ph_amp_map_v1' ) then
#continue
#endif # ! -e
#awk '{print $1,$2,$5}' $ev'_center_ph_amp_map_v1' > $ev'_am.txt_v1'
#mv $ev'_center_ph_amp_map_v1' $ev'.ph.txt_v1'
#C_plot_travel $ev'_am.txt_v1' ~/region_TA 0.1
#C_plot_travel $ev'.ph.txt_v1' ~/region_TA 0.1
#/home/tianye/code/Programs/Earthquake/correct_travel_time_curvature_v1_270 $ev $per 50 0.02 0.05
#rm -f $ev'.ph.txt_v1.HD' $ev'_am.txt_v1.HD'
if( ! -e $ev'.ph.txt_v2' )continue
#awk '{print $1,$2,$5}' $ev'.ph.txt_v2' > $ev'_am.txt_v2'
echo $ev'.ph.txt_v2' $ev'_am.txt_v2' $lon $lat -999. -999. >> ev_files.lst
end

#set para_f='/home/tianye/code/Programs/Attenuation/Alpha_from_ampg_btag_tvtl_cap/param/param_'$per'_'$argv[5]'.txt'
set para_f=$argv[6]
/home/tianye/code/Programs/Attenuation/Alpha_from_ampg_btag_tvtl_cap/alpha_from_ampg_btag_tvtl ev_files.lst $argv[3] $mod_pa'/beta/beta_'$per $para_f

set sm_hdis=`awk '{printf "%.1f",$1}' $para_f`
set p_coef=`awk '{printf "%.1f",$3}' $para_f`
set s_coef=`awk '{printf "%.1f",$4}' $para_f`
set ts=`awk '{printf "%.2f",$12}' $para_f`
set bs=`awk '{printf "%.0f",$13}' $para_f`

#set dir_name = 'test_alp_sm'$sm_hdis'_ts'$ts'_pc'$p_coef'_sc'$s_coef'_'$argv[5]
#if( -e $dir_name ) rm -rf $dir_name
#mv 'test_alp_sm_'$sm_hdis'_'$p_coef'_'$s_coef $dir_name
echo "Smoothing & Ploting..."
set dir_name = 'alpha_sm'$sm_hdis'_ts'$ts'_bs'$bs'_pc'$p_coef'_sc'$s_coef'_'$argv[5]
cd $dir_name
foreach file (alpha_crd_map alpha_raw_map alpha_prp_crd_map alpha_sit_crd_map alpha_site_map alpha_prp_map)
awk '$4>0 && $5>40' $file > $file'_1'
C_plot_travel_nn $file'_1' ~/region_TA 0.1
/home/tianye/code/Programs/Smoothing/Gauss_Smoothing_map_dynamic $file'_1' 70
C_plot_travel_nn $file'_1_smd_70' ~/region_TA 0.1
/home/tianye/code/Programs/Smoothing/Gauss_Smoothing_map_dynamic $file'_1' 100
C_plot_travel_nn $file'_1_smd_100' ~/region_TA 0.1
end
set file='vel_map'
awk '$4>45' $file > $file'_1'
C_plot_travel_nn $file'_1' ~/region_TA 0.1
C_plot_input_region_res_csta $mod_pa'/ph_C/ph_C_'$per'.cpt' $file'_1.HD' ~/region_TA 0.1 0_0
/home/tianye/code/Script/GMT/PLOT/Plot_raw_crd_prp_sit_HGr_HGc.csh $per
/home/tianye/code/Script/GMT/PLOT/Plot_scrd_pcrd_prp_sit_CORp_CORs.csh $per 0
/home/tianye/code/Script/GMT/PLOT/Plot_scrd_pcrd_prp_sit_CORp_CORs.csh $per 1
cd ../..
end

