                   subroutine angles2tensorRad(strike_d,dip_d,slip_d,m)
C       input:  3 angles in degrees
C       output: 6 components of the moment tensor in the order xx,yy,zz,xy,xz,yz
C       Mxx, Mxy, Myy, Mxz, Myz, Mzz
        real*4 m(6)
        data drad/0.0174533/
        strike=strike_d*drad
        dip=dip_d*drad
        slip=slip_d*drad
      m(1) = -sin(dip)*cos(slip)*sin(2.*strike)-sin(2.*dip)*sin(slip)*(sin(strike))**2
      m(2)= sin(dip)*cos(slip)*sin(2.*strike)-sin(2.*dip)*sin(slip)*(cos(strike))**2
      m(3) = sin(2.*dip)*sin(slip)
      m(4)= sin(dip)*cos(slip)*cos(2.*strike)+0.5*sin(2.*dip)*sin(slip)*sin(2.*strike)
      m(5)= -cos(dip)*cos(slip)*cos(strike)-cos(2.*dip)*sin(slip)*sin(strike)
      m(6) = -cos(dip)*cos(slip)*sin(strike)+cos(2.*dip)*sin(slip)*cos(strike)
        return
        end
