foreach st(NZ08)
../CODE/GET_am_ph/whiten_phamp $st"".BDH.SAC 1000 500 5 4
../CODE/GET_am_ph/whiten_phamp $st"".BHZ.SAC 1000 500 5 4
../CODE/lf_get_transfer_function wt$st"".BDH.SAC.am wt$st"".BDH.SAC.ph wt$st"".BHZ.SAC.am wt$st"".BHZ.SAC.ph > TF_$st"".txt
end
