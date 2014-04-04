            subroutine source(sigR,sigL,c1,s1,wvn,v,dvdz,br,bl)
C--to calculate spectral source function bl(6),br(6)
C----------INPUT ARGUMENTS-----------------------------------------
C----code is R or L------------------------------------------------
C----wvn is a wavenumber; c1,s1 are cos,sin of phi (azimuth); 
C----v(3) is an eigenfunction, dvdz is its depth derivative dvdz(3)
C-v1-horizontal,v2-vertical (for Rayleigh),v3 - Love;
C----w-is a frequency in radians;
C---------OUTPUT ARGUMENTS-br (for Rayleigh),bl for Love ----------
C--b(6): b11=b(1);b22=b(2);b33=b(3);b12=b(4);b13=b(5);b23=b(6)
C         xx      yy       zz       xy       xz       yz
C----------------------------------------------------------------------
	    character*1 sigR,sigL
	    real*4 v(3),dvdz(3),wvn(2)
	    complex*8 bl(6),br(6)
C---------step is a spectrum of step-function-or 1 -----------------------
	    do i=1,6
	    bl(i)=(0.0,0.0)
	    br(i)=(0.0,0.0)
	    end do
	    s2=2.*s1*c1
	    csq=c1*c1
	    ssq=s1*s1
	    c2=csq-ssq
	    if(sigR.eq.'+') then   
C---------Rayleigh wave------------
	    aa=-wvn(1)*v(1)
C--------xx------------------------
	    br(1)=cmplx(csq*aa,0.0)
C--------yy------------------------
	    br(2)=cmplx(ssq*aa,0.0)
C--------zz------------------------
	    br(3)=cmplx(dvdz(2),0.0)
C--------xy------------------------
	    br(4)=cmplx(s2*aa,0.0)
	    aa=wvn(1)*v(2)+dvdz(1)
C--------xz------------------------
	    br(5)=cmplx(0.0,aa*c1)
C--------yz------------------------
	    br(6)=cmplx(0.0,aa*s1)
         	end if  
	    if(sigL.eq.'+') then
C--------Love wave-----------------
	    aa=0.5*wvn(2)*v(3)
C--------xx------------------------
	    bl(1)=cmplx(-aa*s2,0.0)
C--------yy------------------------
	    bl(2)=-bl(1)
C--------zz------------------------
            bl(3)=(0.0,0.0)
C--------xy------------------------
	    bl(4)=cmplx(2.*aa*c2,0.0)
C--------xz------------------------
	    aa=dvdz(3)
	    bl(5)=cmplx(0.0,2.*aa*s1)
C--------yz------------------------
	    bl(6)=cmplx(0.0,-2.*aa*c1)
			    end if
			    return
		  end
