c ==========================================================
c Function filter4. Broadband filreting.
c ==========================================================
c Parameters for filter4 function:
c Input parameters:
c dt      - sampling rate of seismogram in seconds, (double)
c n       - number of input samples, (int)
c seis_in - input array length of n, (float)
c Output parameters:
c seis_out - output array length of n, (float)
c ==========================================================

      subroutine fft(dt,n,seis_in,
     1  seis_outamp,seis_outph,nk,dom)
      implicit none
      include 'fftw3.h'
      integer*4 n
      real*8    dt,dt2
      real*4    seis_in(10480576)
      real*4   seis_outamp(10480576), seis_outph(10480576)
c ---
      integer*4 k,ns,nk,i,n2
      real*8    plan1
      real*8    dom
      double complex czero,s(10480576),sf(10480576)
      common    /core/sf,dt2,ns,n2
c ---
      czero = (0.0d0,0.0d0)

c determin the power of FFT
      n2 = n
      ns = 2**max0(int(dlog(dble(n))/dlog(2.0d0))+1,13)
      dt2 = dt
      dom = 1.0d0/dt/ns

      do k = 1,ns
        s(k) = czero
      enddo

      do k = 1,n
        s(k) = seis_in(k)
      enddo
c make backward FFT for seismogram: s ==> sf
      call dfftw_plan_dft_1d(plan1,ns,s,sf, 
     *                  FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1)
      call dfftw_destroy_plan(plan1)
c kill half spectra and correct ends
      nk = ns/2+1
      do k = nk+1,ns
        sf(k) = czero
      enddo
      sf(1) = sf(1)/2.0d0
      sf(nk) = dcmplx(dreal(sf(n)),0.0d0)
c  
       do k = 1,nk
          seis_outamp(k) = 0.0
          seis_outph(k)  = 0.0
       enddo


c      sf(1) = sf(1)/2.0d0
c      sf(nk) = dcmplx(dreal(sf(n)),0.0d0)
      
      do i = 1,nk
         seis_outamp(i)= real(dsqrt(dreal(sf(i))**2 +
     1        dimag(sf(i))**2))
         seis_outph(i) = real(datan2(dimag(sf(i)),dreal(sf(i))))
      enddo
      return
      end

