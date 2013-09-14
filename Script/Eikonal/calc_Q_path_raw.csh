#!/bin/csh
if ( $#argv != 1)then
echo "USAGE: get_event [SAC_path]"
exit
endif
cd $argv[1]

#foreach per ( 10 15 20 25 30 40 50 60 70 80 )
foreach per (25)

@ dis = 4 * $per
cd $per'sec_8snr_'$dis'dis'
rm -f ampf.lst
foreach evinfo (`awk '{print $1"@"$2"@"$3}' ../event.loc`)
set ev=`echo $evinfo | cut -d@ -f1`
if( ! -e $ev'.ph.txt' ) continue
set lon=`echo $evinfo | cut -d@ -f2`
set lat=`echo $evinfo | cut -d@ -f3`
set dist=`~/code/Programs/DIST/get_dist 41 250 $lat $lon d`
if( `echo $dist | awk '{if($1>2000)print 1; else print 0}'` ) continue
if( ! -e $ev'_am.txt' )then
awk '{print $1,$2,$5}' $ev'.ph.txt' > $ev'_am.txt'
endif
echo $ev'_am.txt' >> ampf.lst
end
/home/tianye/code/Programs/Attenuation/alpha_from_amp_ev_lst ampf.lst ../beta/beta_$per ../ph_C/ph_C_$per ../station.lst $per
C_plot_travel_nn 'Western_US_Q_map_'$per'.0sec_raw' ~/region_TA 0.1
C_plot_input_region_res Q_map.cpt 'Western_US_Q_map_'$per'.0sec_raw.HD' ~/region_TA 0.1
awk '{print $1,$2,$4}' 'Western_US_Q_map_'$per'.0sec_raw' > 'Western_US_Q_cfdt_map_'$per'.0sec_raw'
TXT2CPT 'Western_US_Q_cfdt_map_'$per'.0sec_raw'
C_plot_travel_nn 'Western_US_Q_cfdt_map_'$per'.0sec_raw' ~/region_TA 0.1
C_plot_input_region_res 'Western_US_Q_cfdt_map_'$per'.0sec_raw.cpt' 'Western_US_Q_cfdt_map_'$per'.0sec_raw.HD' ~/region_TA 0.1
cd ..
end
