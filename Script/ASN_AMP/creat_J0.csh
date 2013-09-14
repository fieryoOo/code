#foreach a (0.0003 0.0006 0.0009 0.001 0.0015 0.002 0.003)
foreach a (0.0008)
set A0=1.33371e-5
#set a=0.003
set w=`pwd | awk '{print 2*3.1415927/10.889532}'`
set c=3250
awk -v A0=$A0 -v a=$a -v c=$c -v w=$w '{r=$1*c/w; print r,$2*(2*A0**2/a)*exp(-a*r)}' /home/tianye/code/Programs/head/J0.txt > Jtemp_$a'.txt'
end
