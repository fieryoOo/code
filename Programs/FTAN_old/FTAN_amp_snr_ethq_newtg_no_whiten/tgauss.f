c taper phase matched signal
c function [ss] = tgauss(fsnr,gt0,t0,dw,dt,n,seis);
      subroutine tgauss(fsnr,gt0,t0,dw,dt,n,fmatch,seis,
     *                  ss)
      implicit none
      integer*4 n, i,ism,nndl,nndr
      double complex czero,seis(n),ss(n)
      real*8    smax(32768)
c ---
      real*8    pi, dt, gt0, t0, dw, fsnr,sm,smw,fmatch
      integer*4 nc,ii,nnl,nnnl,nnr,nnnr,nleft,nright,left(100),right(100)
      real*8    dzer,dl,dr,vleft(100),vright(100)
      integer*4 nleft_m,nright_m,left_m(100),right_m(100)
      real*8    vleft_m(100),vright_m(100)

      dw = dw
      ism = 1
      dzer = 0.0d0
      czero = (0.0d0,0.0d0)
      pi = datan(1.0d0)*4.0d0
      nc = nint(gt0/dt)+1
c find global max, sm, and index ism
      sm = 0.0d0
      do i = 1,n
         smw = cdabs(seis(i))
         if(smw.ge.sm) then
             sm = smw
             ism = i
         endif
         smax(i) = smw
         ss(i) = seis(i)
      enddo
      write(*,*) 'Distance between maximas=',gt0-(ism-1)*dt-t0,' in sec,',
     * ' Spectra point= ',ism
c find some local minima,# < 100 from left and right side of central max ism
c left side 
      nleft = 0
      do i = ism-1,2,-1     
          dl = smax(i)-smax(i-1)
          dr = smax(i+1)-smax(i)
          if((dl.lt.dzer.and.dr.ge.dzer).or.(dl.le.dzer.and.dr.gt.dzer)) then
              nleft = nleft+1
              left(nleft) = i
              vleft(nleft) = smax(i)
          endif
          if(nleft.ge.100) goto 10
      enddo
   10 continue
c right side
      nright = 0
      do i = ism+1,n-1      
          dl = smax(i)-smax(i-1)
          dr = smax(i+1)-smax(i)
          if((dl.lt.dzer.and.dr.ge.dzer).or.(dl.le.dzer.and.dr.gt.dzer)) then
              nright = nright+1
              right(nright) = i
              vright(nright) = smax(i)
          endif
          if(nright.ge.100) goto 20
      enddo
   20 continue
c find some local maxima,# < 100 from left and right side of central max ism
c left side 
      nleft_m = 0
      do i = ism-1,2,-1
          dl = smax(i)-smax(i-1)
          dr = smax(i+1)-smax(i)
          if((dl.gt.dzer.and.dr.le.dzer).or.(dl.ge.dzer.and.dr.lt.dzer)) then
              nleft_m = nleft_m+1
              left_m(nleft_m) = i
              vleft_m(nleft_m) = smax(i)
          endif
          if(nleft_m.ge.100) goto 110
      enddo
  110 continue
c right side
      nright_m = 0
      do i = ism+1,n-1
          dl = smax(i)-smax(i-1)
          dr = smax(i+1)-smax(i)
          if((dl.gt.dzer.and.dr.le.dzer).or.(dl.ge.dzer.and.dr.lt.dzer)) then
              nright_m = nright_m+1
              right_m(nright_m) = i
              vright_m(nright_m) = smax(i)
          endif
          if(nright_m.ge.100) goto 120
      enddo
  120 continue

c left side, apply cutting
      nnl = 0 
      nnnl = 0
      if(nleft.eq.0) goto 21
       do i = 1,nleft
           if(abs(ism-left(i))*dt.gt.5.0d0) then
                if(vleft(i) .lt. fsnr*sm) then
                    nnl = left(i)
                    nnnl = left_m(i)
                    nndl = nnl-nnnl
                    nndl = min0(nndl,(ism-left(i))/5)
                    nndl = max0(nndl,2)
                    nnl = left(i)-nndl
                    nnnl = left(i)+nndl
                    do ii = 1,nnl
                      ss(ii) = czero
                    enddo
                    do ii = nnl,nnnl
                      ss(ii) = ss(ii)*(1.0d0-dcos(pi/dble(nnnl-nnl)*
     *                (ii-nnl)))/2.0d0
                    enddo
                    goto 21
                endif
           endif
       enddo
c right side, apply cutting
  21  nnr = 0 
      nnnr = 0
      if(nright.eq.0) goto 31
       do i = 1,nright
           if(abs(ism-right(i))*dt.gt.5.0d0) then
                if(vright(i) .lt. fsnr*sm) then
                    nnr = right(i)
                    nnnr = right_m(i)
                    nndr = nnnr-nnr
                    nndr = min0(nndr,(right(i)-ism)/5)
                    nndr = max0(nndr,2)
                    nnr = right(i)-nndr
                    nnnr = right(i)+nndr
                    do ii = nnnr,n
                      ss(ii) = czero
                    enddo
                    do ii = nnr,nnnr
                      ss(ii) = ss(ii)*(1.0d0+dcos(pi/dble(nnnr-nnr)*
     *                (ii-nnr)))/2.0d0
                    enddo
                    goto 31
                endif
           endif
       enddo
   31  continue
      write(*,*) nnl,nnnl,ism,nnr,nnnr
c     stop
      return
      end
