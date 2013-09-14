#!/bin/csh
if ($#argv != 5) then
echo "USAGE: "$0" [amp_sig1] [amp_noise1] [amp_sig2] [amp_noise2] [per]"
exit
endif

set per = $argv[5]
#noise amplitude: 1
set fheader = 'header/ft_2012.FEB.1.J23A.BDH.SAC'
set out1 = '1900.JAN/1900.JAN.1/ft_1900.JAN.1.sig1.BDH.SAC'
set out2 = '1900.JAN/1900.JAN.1/ft_1900.JAN.1.sig2.BDH.SAC'
#produce signal 1
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $fheader $out1 -1 0 $argv[2] 0.
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $out1 $out1 $per 10000 $argv[1] 100.
#/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $out1 $out1 15 10200 3. 150.
#produce signal 2
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $fheader $out2 -1 0 $argv[4] 0.
/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $out2 $out2 $per 10500 $argv[3] 100.
#/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $out2 $out2 15 9700 1. 150.
#/home/tianye/code/Programs/mini_tools/GaussianNoise/ProduceNoise_AddSig $out2 $out2 10 9500 15

~/code/Programs/SEED2COR/Seed2Cor parameters_fake.txt << EOF
y
EOF
cp 1900.JAN/COR/sig1/COR_sig1_sig2.SAC .
