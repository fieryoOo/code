REG=-R0/360/0/0.001
SCA=-JPa5z
psbasemap $REG $SCA -Ba30f10:"Azimuth (degree)":/a0.0002f0.0001:"Amplitude": -X3 -Y3 -P -K > temp1.ps
#makecpt -Z -T0/0.001/0.0002 -D -Cpolar > temp.cpt
awk '{print 0,0,$2,$1,1}' J39A_azi_dis_amp_binavg | psxy $REG $SCA -Svs0.07c/0.5c/0.05c -Ctemp.cpt -O >> temp1.ps

REG=-R0/20/10/25
SCA=-Ja10/15/1/1/1
psbasemap $REG -JX10 -Ba30f10:"Azimuth (degree)":/a2f1:"Amplitude": -X3 -Y3 -P -K > temp2.ps
awk '{print 10,15,$2,$1,1}' J39A_azi_dis_amp_binavg | psxy $REG $SCA -Sv0.07c/0.5c/0.05c -Ctemp.cpt -O >> temp2.ps


