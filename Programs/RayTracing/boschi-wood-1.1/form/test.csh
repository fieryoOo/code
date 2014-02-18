#!/bin/csh
# this script illustrates how everything works
# change incr to get higher-definition images (becomes slower)
set incr = 5
# conversion from generalized harmonic coefficients to xyz files
./gsh2xyz<<EOF
animodel.r100.0
0
$incr
1
40
EOF
./gsh2xyz<<EOF
animodel.r100.2
1
$incr
1
20
EOF
./gsh2xyz<<EOF
animodel.r100.2
2
$incr
1
20
EOF
./gsh2xyz<<EOF
animodel.r100.4
3
$incr
1
20
EOF
./gsh2xyz<<EOF
animodel.r100.4
4
$incr
1
20
EOF
#find fast 2psi or 4psi direction, and azimuthal anisotropy variation
#amplitude, over a grid, from gridded alpha_i files
./xyz2aniso<<EOF
2
animodel.r100.2.eps1.020.xyz
animodel.r100.2.eps2.020.xyz
2
EOF
./xyz2aniso<<EOF
4
animodel.r100.4.eps3.020.xyz
animodel.r100.4.eps4.020.xyz
2
EOF
# use GMT to plot isotropic anomalies
./mapview.csh animodel.r100.0.eps0.040.xyz $incr
# use GMT to plot 2psi fast directions
./plotaniso2.csh animodel.r100.2.eps1.020.xyz.2.x animodel.r100.2.eps1.020.xyz.2.y $incr
# use GMT to plot 4psi fast directions
./plotaniso4.csh animodel.r100.4.eps3.020.xyz.4.x animodel.r100.4.eps3.020.xyz.4.y animodel.r100.4.eps3.020.xyz.4.x2 animodel.r100.4.eps3.020.xyz.4.y2 $incr
# compare the resulting .ps file with those of Boschi and Woodhouse 2006.
# convert back to generalized harmonic coefficients (now only up to L=29)
./xyz2gsh<<EOF
0
animodel.r100.0.eps0.040.xyz
animodel.r100.0.eps0.040.xyz.gsh
0.
EOF
./xyz2gsh<<EOF
2
animodel.r100.2.eps1.020.xyz
animodel.r100.2.eps2.020.xyz
animodel.r100.2.eps.020.xyz.gsh
0.
EOF
./xyz2gsh<<EOF
4
animodel.r100.4.eps3.020.xyz
animodel.r100.4.eps4.020.xyz
animodel.r100.4.eps.020.xyz.gsh
0.
EOF
# now you can compare coefficients in original files animodel.r100.0 animodel.r100.2 and animodel.r100.4 with retrieved coefficients in animodel.r100.*.020.xyz.gsh -- they should be approximately the same (up to L=29 after which we do not retrieve coefficients)
