C######################################################################
        function phaRad (y,x)
c----------------------------------------------------------------------
	pi=3.14159265
	a=sqrt(x**2+y**2)
	if(a.lt.1.e-12) a=1.e-12
	phaRad=asin(y/a)
	if(x.lt.0.)   phaRad=pi-phaRad
C 	if(phaRad.lt.0.) phaRad=phaRad+2.*pi
	RETURN
	END
