c ===========================================================
c ray tracer for epicor program
c ===========================================================
      subroutine tracer(x1,x2,step0,n,t,delta,dbg,ierr)
      implicit none
      integer*4 ierr,n,dbg
      real*8    step0,x1(3),x2(3),t(100)
c ---
      integer*4 iswt,i,nn
      real*8    step,drad,st1,st2,delta,dlen,firot
      real*8    R,xn(3),e(3),vx(4)
c ---
      R = 6371.0d0
      drad = datan(1.0d0)/45.0d0
      call DEL_NORM(x1,x2,xn,delta)
      call rbimod(x1,n,vx,0,ierr)
      if(ierr.ne.0) return
      do i = 1,n
        t(i) = 0.0d0
      enddo
      step = step0*drad
      iswt = 0
      st1 = 0.0d0
      st2 = step
      dlen = 0.0d0
    2 if(st2.ge.delta) then
      st2 = delta
      step = st2-st1
      iswt = 1
      endif
      firot = (st2+st1)/2.0d0
      CALL RTURN(firot,xn,x1,e)
      call NORM(e,e)
      call rbimod(e,nn,vx,dbg,ierr)
      if(ierr.ne.0) return
      if(n.ne.nn) return
      t(1) = t(1)+step/vx(1)
      t(2) = t(2)+step/vx(2)
      t(3) = t(3)+step*vx(3)
      dlen = dlen+step
      if(iswt.eq.1) goto 3
      st1 = st2
      st2 = st2+step 
      goto 2 
    3 do i = 1,nn-1
        t(i) = t(i)*R
      enddo
      call rbimod(x2,nn,vx,0,ierr)
      t(4) = vx(4)
      ierr = 0
      end
