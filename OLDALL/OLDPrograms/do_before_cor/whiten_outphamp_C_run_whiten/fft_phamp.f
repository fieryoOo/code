c ==========================================================
c Function fft_phamp. FFT, output amp and ph files.
c ==========================================================
c Parameters for fft_phamp function:
c Input parameters:
c dt      - sampling rate of seismogram in seconds, (double)
c n       - number of input samples, (int)
c seis_in - input array length of n, (float)
c Output parameters:
c ==========================================================

      subroutine fft_phamp(dt,n,seis_in,
     1  seis_outamp,seis_outph,ns,dom)
      implicit none
      include 'fftw3.h'
      integer*4 n
      real*8    dt
      real*4    seis_in(7500000)
      real*4   seis_outamp(1000000), seis_outph(1000000)
c ---
      integer*4 k,ns,nk,i
      real*8    plan1,plan2
      real*8    dom
      double complex czero,s(7500000),sf(1000000)
c ---
      czero = (0.0d0,0.0d0)



c determin the power of FFT
      ns = 2**((max0(int(dlog(dble(n))/dlog(2.0d0))+1,15))+5.0d0)
      dom = 1.0d0/dt/ns

      do k = 1,ns
        s(k) = czero
      enddo

      do k = 1,n
        s(k) = seis_in(k)
      enddo

c make backward FFT for seismogram: s ==> sf
      call dfftw_plan_dft_1d(plan1,ns,s,sf,
     *                         FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan1)
      call dfftw_destroy_plan(plan1)
c kill half spectra and correct ends

c      nk = ns
      nk = ns/2+1
      do k = nk+1,ns
c      do k = ns/3,ns
        sf(k) = czero
      enddo
 
      sf(1) = sf(1)/2.0d0
      sf(nk) = dcmplx(dreal(sf(n)),0.0d0)

c  
       do k = 1,nk
          seis_outamp(k) = 0.0
          seis_outph(k)  = 0.0
       enddo



       do i = 1,nk
        seis_outamp(i)= real(dsqrt(dreal(sf(i))**2 +
     1                        dimag(sf(i))**2))
        seis_outph(i) = -real(datan2(dimag(sf(i)),dreal(sf(i))))
       enddo


      return
      end

