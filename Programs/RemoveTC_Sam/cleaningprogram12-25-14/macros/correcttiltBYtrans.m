cut off

binoperr npts warning delta warning


sc cd ../J23B

sc cp ../ptcs/J23B.2013.015.dec.ptc ../Platetransfercoeff

setbb azim 28
setbb bh1fn ../J23B/2013.015.00.00.00.0000.7D.J23B..BX1.dec.SAC
setbb bh2fn ../J23B/2013.015.00.00.00.0000.7D.J23B..BX2.dec.SAC
setbb bhzfn ../J23B/2013.015.00.00.00.0000.7D.J23B..BXZ.dec.SAC
setbb outfn ../J23B/2013.015.00.00.00.0000.7D.J23B..BXZ.not.SAC

evaluate to conv2rad 3.14159 / 180.
evaluate to a1 %azim * %conv2rad
evaluate to cos1 cos %a1
evaluate to sin1 sin %a1

# rotate horizontals to tilt direction 
r %bh1fn
rmean
mul %cos1
w temp1c
r %bh2fn
rmean
mul %sin1
addf temp1c
w temph

#predict vertical noise from horizontal record and transferfunction
fft amph
writesp amph temp
sc ./../remveWaterNoise
readsp amph tempf
ifft
rmean
w tempPred

# subtract from vertical
r %bhzfn
subf tempPred
w %outfn

quit

