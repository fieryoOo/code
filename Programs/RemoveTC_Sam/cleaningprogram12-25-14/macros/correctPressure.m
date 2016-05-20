# remove pressure

binoperr npts warning delta warning

sc echo J23B/2013.015.00.00.00.0000.7D.J23B..BXZ.not.SAC

setbb dpgfile J23B/2013.015.00.00.00.0000.7D.J23B..BXH.dec3.SAC
setbb notfile J23B/2013.015.00.00.00.0000.7D.J23B..BXZ.not.SAC
setbb outptrans ../p2zout/transferp2zw.2013.015.J23B
setbb outpcoeff ../ptcp2zs/J23B.2013.015.ptc
setbb outrecord ../J23B/2013.015.00.00.00.0000.7D.J23B..BXZ.nop.SAC


cp ../%dpgfile tempd
cp ../%notfile tempztiltc

sc ./../findP2ZtransferfnsAuto < ../p2z.w.inp
cp tempP2Ztransfer %outptrans
cp tempP2Zcoeff %outpcoeff
cut off
r tempd
rmean
taper
# tapering because very long-period tidal components
# cause discontinuities at ends of record
fft amph
writesp amph temp
sc ./../remveWaterNoiseAuto
readsp amph tempf
ifft
rmean
w tempPred
# enter vertical filename
r tempztiltc
rmean
taper
subf tempPred
w %outrecord

quit

