#!/bin/csh
python << END
t0=15*60+19.723
dis_fvm=2668.5
dis_usin=2910.3
inter_dis=dis_usin-dis_fvm
FVM=[[27,13.9],[27,37.2],[28,3.2],[28,23.2],[28,39.9],[28,57.3],[29,11.3],[29,27.3],[29,39.9],[29,56.6],[30,7.3],[30,16.7]]
USIN=[[28,11.9],[28,39.3],[29,3.9],[29,23.9],[29,42.7],[30,1.3],[30,17.3],[30,30.7],[30,44.1],[30,58],[31,12.7],[31,23.4]]
fvm=[]
usin=[]
for i in range(12):
  fvm.append(FVM[i][0]*60+FVM[i][1]-t0)
  usin.append(USIN[i][0]*60+USIN[i][1]-t0)

fvm_dsp=[[] for i in range(10)]
usin_dsp=[[] for i in range(10)]
ph_dsp=[[] for i in range(10)]
for i in range(1,11):
  fvm_dsp[i-1].append(fvm[i+1]-fvm[i-1])
  fvm_dsp[i-1].append(dis_fvm/fvm[i])
  usin_dsp[i-1].append(usin[i+1]-usin[i-1])
  usin_dsp[i-1].append(dis_usin/usin[i])
  ph_dsp[i-1].append((fvm_dsp[i-1][0]+usin_dsp[i-1][0])/2)
  ph_dsp[i-1].append(inter_dis/(usin[i]-fvm[i]))
  print i,usin[i]-fvm[i]
ff = open("fvm_group",'w')
for i in range(10):
  ff.write("%.2f   %.3f\n" % (fvm_dsp[i][0],fvm_dsp[i][1]))

ff.close()

ff = open("usin_group",'w')
for i in range(10):
  ff.write("%.2f   %.3f\n" % (usin_dsp[i][0],usin_dsp[i][1]))

ff.close()

ff = open("ph_dsp",'w')
for i in range(10):
  ff.write("%.2f   %.3f\n" % (ph_dsp[i][0],ph_dsp[i][1]))

ff.close()
END

gnuplot << END
set xlabel "Period (sec)"
set ylabel "V (km/sec)"
set xrange [1:100]
set yrange [2.5:4.5]
set ytics 2.4,0.4
set mytics 4
set logscale x

plot 'fvm_group' using 1:2 title 'Rayleigh Group Velocity - Station FVM'
set term postscript enhanced color
set output 'fvm_group.ps'
replot
set output
set term x11

plot 'usin_group' using 1:2 title 'Rayleigh Group Velocity - Station USIN'
set term postscript enhanced color
set output 'usin_group.ps'
replot
set output
set term x11

set ylabel "C (km/sec)"
plot 'ph_dsp' using 1:2 with line title 'Rayleigh Phase Velocity - Station FVM_USIN'
set term postscript enhanced color
set output 'ph_dsp.ps'
replot
set output
set term x11
END
