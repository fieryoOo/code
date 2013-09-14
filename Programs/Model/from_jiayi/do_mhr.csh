#/bin/csh
#execute mhr to get the predicted travel time for each path

set fpath = ./pathfile
#set vel = /home/jiayi/Tool/MODEL/SURF_VEL/PHASE_VEL/BIGCHINA/cmb_phv
set vel = /home/jiayi/Tool/MODEL/SURF_VEL/PHASE_VEL/WESTCHINA/cmb2_phv
set perlst = ./perlist
set dir = PREDICTION
if ( 1 == 0 )then
if ( ! -e  $dir ) then
	mkdir $dir
else
	rm $dir
	mkdir $dir
	echo "exist already"
endif

echo $fpath $vel $perlst

/home/jiayi/progs/jy/GLOBAL_DISP/mhr_grvel_predict/mhr_grvel_predict_earth_v3_cv_for_sm $fpath $vel $perlst
mv -f PREDICTION_R $dir
mv -f PREDICTION_L $dir

endif
set cmp = L
if (! -e $dir/evt_sta_pair) then
mkdir $dir/evt_sta_pair
mkdir $dir/evt_sta_pair/$cmp
endif
/home/jiayi/progs/jy/GLOBAL_DISP/get_COR_earth $dir/PREDICTION_$cmp $dir/evt_sta_pair/$cmp >& err_log_cor_pred

