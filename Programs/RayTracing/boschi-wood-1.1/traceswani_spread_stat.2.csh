#!/bin/csh
set iorbit = 1
set fileout = spread
set c0 = 4.08957 # phase velocity of Rayleigh at 100s in PREM
set prd = 100. # period 
./traceswani_spread_stat.2<<EOF
$c0
$prd
animodel.r100.0
animodel.r100.2
animodel.r100.4
$iorbit
$fileout
random_locs.txt
EOF
