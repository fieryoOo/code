#!/bin/csh
#echo on
if ($#argv != 2) then
  echo "USAGE: tool_do_AtoD [region_infile] [sac_path : sac_path/sec_snr_dis/evt.ph.txt]"
  exit 1
endif
#/home/tianye/data_Eikonal/SAC_XR/80sec_10snr_960dis/temp
#/home/tianye/data_Eikonal/SAC_TA/80sec_10snr_960dis/#/home/tianye/data_Eikonal/SAC_XR/80sec_10snr_960dis/temp
set here = `pwd`

cd $argv[2]

set REG = `more $argv[1]`
#echo $REG
set sta_num_cri = 50
set cur_tm_cri = 0.005
set cur_amp_cri = 0.0625 #1/16   
set cla = 30.0
set clo = 110.0
set grid = 0.1
set snr = 10
set distlst = ( 0 15000 )
set here = `pwd`
foreach per ( 10 15 20 )
    foreach i(  1  )
	set j = `echo $i | awk '{print $1+1}'`
	set dis1 = ${distlst[$i]}
	set dis2 = ${distlst[$j]}
	#cd $per"sec_${snr}snr_${dis1}_${dis2}dist_noAMPcor"
	cd AMP_COR_dd150da4_phtxt_Nd10_Na6_all_renewcode/$per"sec_${snr}snr_${dis1}_${dis2}dist"
	echo "=====WORKING ON PER="$per" ================"
##########################################
#rm -f ampwrong_* pha2pi_* phawrong_* *_v1
#rm -f *ph.txt_v1.HD *_am.txt_v1.HD cur_tm_big* cur_amp_big* *.ph.txt_v2* *_am.txt_v2*
############################################
        #rm -f *.ph.txt_v1*  *.ph.txt_v2*
	#foreach event(`ls 20*.ph.txt | cut -d. -f1`)
	#foreach event ( 20090629180350  20090529062014 20090401035458 )
	foreach event ( 20090629180350 )
		#rm -f ampwrong_*$event*  pha2pi_*$event*  phawrong_*event*
		#echo "get v1"
		echo "-------------DO phase"	
		/home/jiayi/progs/jy/eikonal/earthquake/correct_2pi_v1_jy_ph $event $per $sta_num_cri $cla $clo
		/home/jiayi/Script/GMT/C_plot_travel_g0.4 $event".ph.txt_v1" $here/$argv[1]  $grid
		/home/jiayi/progs/jy/eikonal/earthquake/correct_travel_time_curvature_v1_270_jy_ph_g0.4 $event $per $sta_num_cri $cur_tm_cri $grid
		rm -f `wc $event'.ph.txt_v2' | grep "0       0       0" | awk '{print $4}'`
		/home/jiayi/Script/GMT/C_plot_travel_g0.4 $event'.ph.txt_v2' $here/$argv[1] $grid
		/home/jiayi/Script/GMT/C_plot_travel_T0.2_g0.4 $event'.ph.txt_v2' $here/$argv[1] $grid
		awk '{split($0,a," ");printf"%8.1f%8.1f\t",a[1],a[2];getline < "'$event'.ph.txt_v2.HD_0.2";printf" %10f\t%10f\t%10f\n", a[3]-$3,a[3],$3}' $event.ph.txt_v2.HD  > $event.ph.txt_v2.HD_dt
		echo "---------DO amp"
		/home/jiayi/progs/jy/eikonal/earthquake/correct_2pi_v1_jy_am $event $per $sta_num_cri $cla $clo
		#/home/jiayi/progs/jy/eikonal/earthquake/correct_2pi_v1_jy_am_tmp $event $per $sta_num_cri $cla $clo
		rm -f  ${event}_ph.txt_v1.HD ${event}_am.txt_v1.HD  cur_tm_big_${event}_v1.txt cur_amp_big_${event}_v1.txt 
		/home/jiayi/Script/GMT/C_plot_travel_amv1_jy_g0.4  $event".am.txt_v1" $here/$argv[1] $grid
		/home/jiayi/progs/jy/eikonal/earthquake/correct_travel_time_curvature_v1_270_jy_am_g0.4 $event $per $sta_num_cri  $cur_amp_cri $grid
		rm -f `wc $event'.am.txt_v2' | grep "0       0       0" | awk '{print $4}'`
		/home/jiayi/Script/GMT/C_plot_travel_amv2_jy_g0.4 $event'.am.txt_v2' $here/$argv[1] $grid
	

	end #foreach event
	cd $here
    end #foreach i
end #foreach per
