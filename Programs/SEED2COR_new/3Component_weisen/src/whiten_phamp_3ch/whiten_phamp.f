c ==========================================================
c Function filter4. Broadband filreting.
c ==========================================================
c Parameters for filter4 function:
c Input parameters:
c f1,f2   - low corner frequences, f2 > f1, Hz, (double)
c f3,f4   - high corner frequences, f4 > f3, Hz, (double)
c npow    - power of cosine tapering,  (int)
c dt      - sampling rate of seismogram in seconds, (double)
c n       - number of input samples, (int)
c seis_in - input array length of n, (float)
c Output parameters:
c seis_out - output array length of n, (float)
c ==========================================================

      subroutine filter4(f1,f2,f3,f4,npow,dt,n,seis_inZ,
     1  seis_outZ,seis_outampZ,seis_outphZ,
     1  seis_inE,seis_outE,seis_outampE,seis_outphE,
     1  seis_inN,seis_outN,seis_outampN,seis_outphN,
     1  ns,dom)
      implicit none
      include '../inc/fftw3.h'
      integer*4 npow,n
      real*8    f1,f2,f3,f4,dt,temp
      real*4    seis_inZ(400000),seis_outZ(400000)
      real*4   seis_outampZ(400000), seis_outphZ(400000)
      real*4    seis_inE(400000),seis_outE(400000)
      real*4   seis_outampE(400000), seis_outphE(400000)
      real*4    seis_inN(400000),seis_outN(400000)
      real*4   seis_outampN(400000), seis_outphN(400000)
c ---
      integer*4 k,ns,nk,i
      real*8    plan1,plan2
      real*8    dom,dpi
      double complex czero,sZ(400000),sfZ(400000)
      double complex sE(400000),sfE(400000)
      double complex sN(400000),sfN(400000)
c ---
      czero = (0.0d0,0.0d0)
      dpi = datan(1.0d0)*4.0d0


c determin the power of FFT
      ns = 2**max0(int(dlog(dble(n))/dlog(2.0d0))+1,13)
      dom = 1.0d0/dt/ns

      do k = 1,ns
        sZ(k) = czero
        sE(k) = czero
        sN(k) = czero
      enddo

      do k = 1,n
        sZ(k) = seis_inZ(k)
        sE(k) = seis_inE(k)
        sN(k) = seis_inN(k)
      enddo

c make backward FFT for seismogram: s ==> sf
      call dfftw_plan_dft_1d(plan1,ns,sZ,sfZ,
     *                         FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1)
      call dfftw_destroy_plan(plan1)
      call dfftw_plan_dft_1d(plan1,ns,sE,sfE,
     *                         FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1)
      call dfftw_destroy_plan(plan1)
      call dfftw_plan_dft_1d(plan1,ns,sN,sfN,
     *                         FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1)
      call dfftw_destroy_plan(plan1)
c kill half spectra and correct ends

      nk = ns/2+1
      do k = nk+1,ns
        sfZ(k) = czero
        sfE(k) = czero
        sfN(k) = czero
      enddo
 
      sfZ(1) = sfZ(1)/2.0d0
      sfZ(nk) = dcmplx(dreal(sfZ(n)),0.0d0)
      sfE(1) = sfE(1)/2.0d0
      sfE(nk) = dcmplx(dreal(sfE(n)),0.0d0)
      sfN(1) = sfN(1)/2.0d0
      sfN(nk) = dcmplx(dreal(sfN(n)),0.0d0)
c  
       do k = 1,nk
          seis_outampZ(k) = 0.0
          seis_outphZ(k)  = 0.0
          seis_outampE(k) = 0.0
          seis_outphE(k)  = 0.0
          seis_outampN(k) = 0.0
          seis_outphN(k)  = 0.0
       enddo


c=============================================================
c     do smoothing on sf equivalent to do " smooth mean h 20" in SAC

      call smooth(f1,f2,f3,f4,dom,nk,sfZ,sfE,sfN,20)
C=============================================================
c===============================================================
c   make tapering
      call flt4(f1,f2,f3,f4,dom,nk,npow,sfZ,sfE,sfN)

c       do i = 1,nk
c        seis_outamp(i)= real(dsqrt(dreal(sf(i))**2 +
c     1                        dimag(sf(i))**2))
c        seis_outph(i) = real(datan2(dimag(sf(i)),dreal(sf(i))))
c       enddo


c make forward FFT for seismogram: sf ==> s
      call dfftw_plan_dft_1d(plan2,ns,sfZ,sZ,
     *                         FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan2)
      call dfftw_destroy_plan(plan2)
      call dfftw_plan_dft_1d(plan2,ns,sfE,sE,
     *                         FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan2)
      call dfftw_destroy_plan(plan2)
      call dfftw_plan_dft_1d(plan2,ns,sfN,sN,
     *                         FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan2)
      call dfftw_destroy_plan(plan2)
c forming final result
      do k = 1,n
         seis_outZ(k) = 2.0*real(dreal(sZ(k)))/ns
         seis_outE(k) = 2.0*real(dreal(sE(k)))/ns
         seis_outN(k) = 2.0*real(dreal(sN(k)))/ns
      enddo
      do k= 1,500
         temp=(1-dcos((k-1.0)/500.0*dpi))/2.0d0 
         seis_outZ(k)=seis_outZ(k)*temp
         seis_outZ(n+1-k)=seis_outZ(n+1-k)*temp
         seis_outE(k)=seis_outE(k)*temp
         seis_outE(n+1-k)=seis_outE(n+1-k)*temp
         seis_outN(k)=seis_outN(k)*temp
         seis_outN(n+1-k)=seis_outN(n+1-k)*temp
         
      enddo
      do k = n+1,ns
         sZ(k) = czero
         sE(k) = czero
         sN(k) = czero
      enddo
      do k = 1,n
        sZ(k) = seis_outZ(k)
        sE(k) = seis_outE(k)
        sN(k) = seis_outN(k)
      enddo
c       make backward FFT for seismogram: s ==> sf
      call dfftw_plan_dft_1d(plan1,ns,sZ,sfZ,
     *                         FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1)
      call dfftw_destroy_plan(plan1)
      call dfftw_plan_dft_1d(plan1,ns,sE,sfE,
     *                         FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1)
      call dfftw_destroy_plan(plan1)
      call dfftw_plan_dft_1d(plan1,ns,sN,sfN,
     *                         FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1)
      call dfftw_destroy_plan(plan1)
c      do k = nk+1,ns
c         sf(k) = czero
c      enddo
      
      sfZ(1) = sfZ(1)/2.0d0
      sfZ(nk) = dcmplx(dreal(sfZ(n)),0.0d0)
      sfE(1) = sfE(1)/2.0d0
      sfE(nk) = dcmplx(dreal(sfE(n)),0.0d0)
      sfN(1) = sfN(1)/2.0d0
      sfN(nk) = dcmplx(dreal(sfN(n)),0.0d0)
c     
c      do k = 1,nk
c         seis_outamp(k) = 0.0
c         seis_outph(k)  = 0.0
c      enddo
      do i = 1,nk
         seis_outampZ(i)= real(dsqrt(dreal(sfZ(i))**2 +
     1        dimag(sfZ(i))**2))
         seis_outphZ(i) = real(datan2(dimag(sfZ(i)),dreal(sfZ(i))))
         seis_outampE(i)= real(dsqrt(dreal(sfE(i))**2 +
     1        dimag(sfE(i))**2))
         seis_outphE(i) = real(datan2(dimag(sfE(i)),dreal(sfE(i))))
         seis_outampN(i)= real(dsqrt(dreal(sfN(i))**2 +
     1        dimag(sfN(i))**2))
         seis_outphN(i) = real(datan2(dimag(sfN(i)),dreal(sfN(i))))
      enddo
      return
      end



c===============================================================
c Tapering subroutine itself
c===============================================================
      subroutine flt4(f1,f2,f3,f4,dom,nk,npow,sfZ,sfE,sfN)
      real*8    f1,f2,f3,f4,dom
      integer*4 nk,npow
      double complex sfZ(400000),sfE(400000),sfN(400000)
      real*8    d1,d2,f,dpi,ss,s(400000)
      integer*4 i,j
c ---
      dpi = datan(1.0d0)*4.0d0
      do i = 1,nk
         s(i) = 0.0d0
      enddo
      do i = 1,nk
        f = (i-1)*dom
        if(f.le.f1) then
          goto 1
        else if(f.le.f2) then
          d1 = dpi/(f2-f1)
          ss = 1.0d0
          do j = 1,npow
            ss = ss*(1-dcos(d1*(f1-f)))/2.0d0
          enddo
          s(i) = ss
        else if(f.le.f3) then
           s(i) = 1.0d0
        else if(f.le.f4) then
          d2 = dpi/(f4-f3)
          ss = 1.0d0
          do j = 1,npow
            ss = ss*(1+dcos(d2*(f3-f)))/2.0d0
          enddo
          s(i) = ss
        endif
  1     continue
      enddo
      do i = 1,nk
        sfZ(i) = sfZ(i)*s(i)
        sfE(i) = sfE(i)*s(i)
        sfN(i) = sfN(i)*s(i)
      enddo
      return
      end


c===============================================================
c  smoothing routine      call smooth(f1,f2,f3,f4,dom,nk,sf,20)
c=s==============================================================
      subroutine smooth(f1,f2,f3,f4,dom,nk,sfZ,sfE,sfN,number)
      real*8    f1,f2,f3,f4
      integer*4 number,nk
      double complex sfZ(400000),sfE(400000),sfN(400000)
      real*8    sorigZ(400000), sout(400000),dom
      real*8    sorigE(400000) 
      real*8    sorigN(400000)
      real*8   f,sum, avg
c ---
        do i = 1,nk
         sorigZ(i) = dsqrt(dreal(sfZ(i))**2+dimag(sfZ(i))**2)
         sorigE(i) = dsqrt(dreal(sfE(i))**2+dimag(sfE(i))**2)
         sorigN(i) = dsqrt(dreal(sfN(i))**2+dimag(sfN(i))**2)
        enddo
     
        do i = 1,nk

        f = (i-1)*dom

        if( f .ge. f1 .and. f .le. f4 ) then
            sum = 0. 
          do jk = -number,number
             ijk = i+jk
             sum = sum + sorigZ(ijk)
             sum = sum + sorigE(ijk)
             sum = sum + sorigN(ijk)
          enddo
            sout(i) = sum/(2.*number+1.)/3
        else
            sout(i) = (sorigZ(i)+sorigE(i)+sorigN(i))/3
        endif

       enddo


       do i = 1,nk
         f = (i-1)*dom
       if( f .ge. f1 .and. f .le. f4 ) then
          sout(i) = 1.0d0/sout(i)
       else
          sout(i) = 0.0d0
       endif
          
       enddo



        do i = 1,nk
           sfZ(i) = sfZ(i)*sout(i)
           sfE(i) = sfE(i)*sout(i)
           sfN(i) = sfN(i)*sout(i)
        enddo

       return
 
       end


