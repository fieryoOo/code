c******************************************************
c Rotation around normal vector
c******************************************************
      SUBROUTINE RTURN(fi,a,x,x1)
      IMPLICIT NONE
      real*8    a(3),x(3),x1(3),r(3),r1(3),co,si,ax,fi
      integer*4 i
      co=DCOS(fi)
      si=DSIN(fi)
      ax=a(1)*x(1)+a(2)*x(2)+a(3)*x(3)
      do 10 i=1,3
   10 r(i)=x(i)-a(i)*ax
      CALL VECT(r,a,r1)
      do 20 i=1,3
   20 x1(i)=r(i)*co-r1(i)*si +a(i)*ax
      return
      end
c******************************************************
c Normal vector and distance evaluation
c******************************************************
      SUBROUTINE DEL_NORM(x1,x2,xn,delta)
      IMPLICIT NONE
      real*8 x1(3),x2(3),xn(3),delta
      CALL NORM(x1,x1)
      CALL NORM(x2,x2)
      CALL VECT(x1,x2,xn)
      CALL NORM(xn,xn)
      delta=DACOS(x1(1)*x2(1)+x1(2)*x2(2)+x1(3)*x2(3))
      return
      end
c******************************************************
c Azimuth evaluation.
c******************************************************
      SUBROUTINE SAZI(e,xn,az)
      implicit none
      real*8 e(3),xn(3),azn(3),azpole(3),az,dpi
      real*8 d,a(3),b(3),c(3)
      data azpole /0.0d0,0.0d0,1.0d0/
      dpi=DATAN(1.d0)*4.d0
      CALL NVECT(azpole,e,c,azn,d)
      CALL NVECT(e,azn,c,b,d)
      CALL NVECT(xn,e,c,a,d)
      CALL VECT(b,a,c)
      d=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      if(DABS(d).gt.1.0d0) d= DSIGN(1.0d0,d)
      az=DACOS(d)
      if((e(1)*c(1)+e(2)*c(2)+e(3)*c(3)).gt.0.0d0) az=2.0d0*dpi-az
      return
      end
c******************************************************
c Vector product: c=[a x b]
c******************************************************
      SUBROUTINE VECT(a,b,c)
      implicit none
      real*8    a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end
c******************************************************
c Vector normalization: b=a/||a||
c******************************************************
      SUBROUTINE NORM(a,b)
      implicit none
      real*8    a(3),b(3),dnorm
      integer*4 i
      dnorm=DSQRT(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      do 10 i=1,3
   10 b(i)=a(i)/dnorm
      return
      end
c******************************************************
c Vector product: c=[a x b], d=||c||,cn=c/||c||
c******************************************************
      SUBROUTINE NVECT(a,b,c,cn,d)
      implicit none
      integer*4 i
      real*8    a(3),b(3),c(3),cn(3),d
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      d=DSQRT(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
      if(d.lt.1.0d-15) then
      if(a(1).ne.0.0d0) then
      cn(1)=-a(2)
      cn(2)=a(1)
      cn(3)=a(3)
      else
      cn(1)=a(1)
      cn(2)=-a(3)
      cn(3)=a(2)
      endif
      d=DSQRT(cn(1)*cn(1)+cn(2)*cn(2)+cn(3)*cn(3))
      do 10 i=1,3
   10 cn(i)=cn(i)/d
      d=0.0d0
      return
      endif
      do 20 i=1,3
   20 cn(i)=c(i)/d
      return
      end
c******************************************************
c  c = [a x b], cn = c / ||c||, d = ||c||
c******************************************************
      subroutine svect(a,b,c,cn,d)
      implicit none
      real*4    d,a(3),b(3),c(3),cn(3)
      integer*4 i
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      d=sqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
      if(d.eq.0.0d0) then
      if(a(1).ne.0.0d0) then
      cn(1)=-a(2)
      cn(2)=a(1)
      cn(3)=a(3)
      else
      cn(1)=a(1)
      cn(2)=-a(3)
      cn(3)=a(2)
      endif
      d=sqrt(cn(1)*cn(1)+cn(2)*cn(2)+cn(3)*cn(3))
      do 10 i=1,3
   10 cn(i)=cn(i)/d
      d=0.0d0
      return
      endif
      do 20 i=1,3
   20 cn(i)=c(i)/d
      return
      end
