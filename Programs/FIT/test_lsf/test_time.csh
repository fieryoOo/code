set i = 2
while ( $i < 20 )
set N = `echo $i | awk '{print 2**$1}'`
pwd | awk -v N=$N '{for(i=0;i<N;i++)print i,3*i+10+(rand()-0.5)*10.}' > temp.dat
set tmp = `least_squares_line temp.dat 0`
set result = `echo $tmp | awk '{print $1,$2}'`
set tm = `echo $tmp | awk '{print $3}'`
echo $N $tm" ( "fit results: $result" )"
@ i ++
end
