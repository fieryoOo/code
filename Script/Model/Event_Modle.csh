#/bin/csh
#execute mhr to get the predicted travel time for each path
if ( $#argv != 4)then
echo "USAGE: Event_Model.csh [station.lst] [event.loc] [phvel_map_name] [out_dir]"
exit
endif

cp $argv[1] /home/tianye/code/Programs/Model/station.lst
cp $argv[2] /home/tianye/code/Programs/Model/event.loc
cd /home/tianye/code/Programs/Model/
python /home/tianye/code/Programs/Model/get_input_hmr_sta.py
mv pathfile /home/tianye/Model/PHASE_map/
cd /home/tianye/Model/PHASE_map/
awk '$3!=$4' pathfile > temp1
mv temp1 pathfile
set fpath = ./pathfile
#set vel = /home/jiayi/Tool/MODEL/SURF_VEL/PHASE_VEL/WESTCHINA/cmb2_phv
#set vel = 'USglobal'
set vel = $argv[3]
set perlst = ./per.lst
set dir = $argv[4]

echo $fpath $vel $perlst

mkdir -p $dir
/home/tianye/code/Programs/Model/mhr_grvel_predict_earth_v2 $fpath $vel $perlst
mv -f PREDICTION_R $dir
mv -f PREDICTION_L $dir
cd $dir
/home/tianye/Model/PHASE_map/extract_PREDICTION PREDICTION_R
exit

set cmp = L
if (! -e $dir/evt_sta_pair) then
mkdir $dir/evt_sta_pair
mkdir $dir/evt_sta_pair/$cmp
endif
/home/jiayi/progs/jy/GLOBAL_DISP/get_COR_earth $dir/PREDICTION_$cmp $dir/evt_sta_pair/$cmp >& err_log_cor_pred

