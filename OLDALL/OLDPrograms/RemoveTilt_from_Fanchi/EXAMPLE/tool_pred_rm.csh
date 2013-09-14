../CODE/PRED_OBS_from_DPG/pred_from_transfer_function NZ08.BDH.SAC ../Quiet/TF_NZ08.txt NZ08.PRED.SAC 
sac<<END
r NZ08.BHZ.SAC 
subf  NZ08.PRED.SAC 
w NZ08.CORRECTED.SAC
quit
END
