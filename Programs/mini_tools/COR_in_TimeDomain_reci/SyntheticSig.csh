#noise amplitude: 1
set fheader = '2012.FEB/2012.FEB.1/ft_2012.FEB.1.J23A.BDH.SAC'
set out1 = '1900.JAN/1900.JAN.1/ft_1900.JAN.1.sig1.BDH.SAC'
set out2 = '1900.JAN/1900.JAN.1/ft_1900.JAN.1.sig2.BDH.SAC'
#produce signal 1
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $fheader $out1 -1 0 0
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $out1 $out1 10 10000 10
#produce signal 2
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $fheader $out2 -1 0 0
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $out2 $out2 10 10100 5
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $out2 $out2 10 9800 15
