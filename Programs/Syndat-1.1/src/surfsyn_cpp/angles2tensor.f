                   subroutine angles2tensor(the_d,dip_d,slip_d,m)
C       input:  3 angles in degrees
C       output: 6 components of the moment tensor in the order xx,yy,zz,xy,xz,yz
C       Mxx, Myy, Mzz, Mxy, Mxz, Myz
        real*4 m(6)
        data drad/0.0174533/
        the=the_d*drad
        dip=dip_d*drad
        slip=slip_d*drad
      m(1) = -sin(dip)*cos(slip)*sin(2.*the)-sin(2.*dip)*sin(slip)*(sin(the))**2
      m(2)= sin(dip)*cos(slip)*sin(2.*the)-sin(2.*dip)*sin(slip)*(cos(the))**2
      m(3) = sin(2.*dip)*sin(slip)
      m(4)= sin(dip)*cos(slip)*cos(2.*the)+0.5*sin(2.*dip)*sin(slip)*sin(2.*the)
      m(5)= -cos(dip)*cos(slip)*cos(the)-cos(2.*dip)*sin(slip)*sin(the)
      m(6) = -cos(dip)*cos(slip)*sin(the)+cos(2.*dip)*sin(slip)*cos(the)
C      do i=1,6
C      m(i)=m(i)*m0
C      enddo
        return
        end
C     character*10 sym
C     real*4 m(6),m0
C     call Getarg(1,sym)
C     read(sym,*) the_d
C     call Getarg(2,sym)
C     read(sym,*) dip_d
C     call Getarg(3,sym)
C     read(sym,*) slip_d
C     call Getarg(4,sym)
C     read(sym,*) m0     
C     call angles2tensor(the_d,dip_d,slip_d,m,m0)
C     PRINT*,m
C     end
