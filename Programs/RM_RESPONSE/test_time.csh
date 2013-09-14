### warmup
set i = 1
while( $i < 5000 )
echo AAA >& /dev/null
@ i ++
end
### initialize
set SACf = 2012.FEB/2012.FEB.11/2012.FEB.11.J23A.BHZ.SAC
set RESf = 2012.FEB/2012.FEB.11/RESP.7D.J23A..BHZ
set i = 0
### run RemoveRESP 50 times ( result: 20.65 sec )
if( 0 ) then
time
while( $i < 50 )
RemoveRESP $SACf $RESf tmp.sac
@ i ++
echo $i
end
time
endif
### run SAC 50 times ( result: 107.78 sec )
if( 1 ) then
time
while( $i < 50 )
sac << END
r $SACf
rmean
rtr
transfer from EVALRESP FNAME $RESf to vel freqlimits 0.01125 0.0125 1.0 1.1
w tmp.sac
quit
END
@ i ++
end
time
endif
